library(data.table)
library(sp)
library(spBayes)
library(coda)

# ==============================================================================
# [설정] 경로 및 공통 파라미터
# ==============================================================================

OUT_RD_DIR       <- "C:/Users/USER/Desktop/thesis/output1/sgg_rd"
OUT_BAYES_DIR    <- file.path(OUT_RD_DIR, "bayes_result_sensitivity")
if (!dir.exists(OUT_BAYES_DIR)) dir.create(OUT_BAYES_DIR, recursive = TRUE)

years <- 2015:2024

# MCMC 설정 (방법 설명에 맞게)
n.samples <- 10000
burn.in   <- 7500

# 초기값 및 튜닝 (두 시나리오에서 공통 사용)
starting <- list("phi" = 1, "sigma.sq" = 5, "tau.sq" = 1)
tuning   <- list("phi" = 0.1, "sigma.sq" = 0.1, "tau.sq" = 0.1)

# ==============================================================================
# [함수] 한 번의 spLM + spRecover + SNR 계산을 수행하는 헬퍼
# ==============================================================================

run_spLM_with_priors <- function(bayesDF, coords, priors, scenario_name,
                                 n.samples, burn.in,
                                 starting, tuning, cov.model = "spherical",
                                 n.report = 2000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  cat("    · 시나리오:", scenario_name, " / n =", nrow(bayesDF), "/n")
  
  # 모델 피팅
  m.1 <- spLM(
    month_wt_srd ~ 1,
    coords    = coords,
    data      = bayesDF,
    starting  = starting,
    tuning    = tuning,
    priors    = priors,
    cov.model = cov.model,
    n.samples = n.samples,
    verbose   = FALSE,
    n.report  = n.report
  )
  
  # theta (분산 관련 파라미터들) posterior sample
  theta.samples <- m.1$p.theta.samples
  if (burn.in >= 1) {
    theta.samples <- theta.samples[burn.in:n.samples, , drop = FALSE]
  }
  
  # sigma.sq, tau.sq 추출 후 SNR 계산
  sigma.sq.post <- theta.samples[, "sigma.sq"]
  tau.sq.post   <- theta.samples[, "tau.sq"]
  snr.post      <- sigma.sq.post / (sigma.sq.post + tau.sq.post)
  
  snr.mean <- mean(snr.post, na.rm = TRUE)
  snr.sd   <- sd(snr.post, na.rm = TRUE)
  
  # 공간 효과 복원
  m.1.rec <- spRecover(m.1, start = burn.in, verbose = FALSE)
  
  beta.samples <- m.1.rec$p.beta.recover.samples
  w.samples    <- m.1.rec$p.w.recover.samples
  
  beta.mean <- mean(beta.samples)
  w.mean    <- apply(w.samples, 1, mean)
  w.sd      <- apply(w.samples, 1, sd)
  
  # 결과 리스트로 반환
  list(
    beta.mean  = beta.mean,
    w.mean     = w.mean,
    w.sd       = w.sd,
    snr.mean   = snr.mean,
    snr.sd     = snr.sd,
    theta.post = theta.samples  # 필요하면 나중에 따로 활용 가능
  )
}

# ==============================================================================
# [루프] 연도별로 Flat vs Informative prior 비교
# ==============================================================================

# 요약 결과 저장용
summary_list <- list()

for (YR in years) {
  cat("==========================================================/n")
  cat("▶", YR, "년 베이지안 공간 민감도 분석 시작/n")
  cat("==========================================================/n")
  
  tryCatch({
    
    # 1. 데이터 로드
    file_path <- file.path(OUT_RD_DIR, sprintf("%d_sgg_rd_monthwt.csv", YR))
    if (!file.exists(file_path)) {
      cat("  [Skip]", YR, "년 입력 파일이 없습니다./n")
      next
    }
    
    baz <- fread(file_path, encoding = "UTF-8")
    
    # 2. 전처리
    bayesDF <- baz[, .(SGG_CD, month_wt_srd, lon, lat, pop)]
    bayesDF <- bayesDF[
      !is.na(month_wt_srd) & is.finite(month_wt_srd) &
        !is.na(lon) & !is.na(lat) & !is.na(pop)
    ]
    
    # 인구 필터링 (기존 로직 유지)
    bayesDF <- bayesDF[pop > 1000]
    
    if (nrow(bayesDF) < 10) {
      cat("  [Skip] 유효 데이터 부족 (n <", nrow(bayesDF), ")./n")
      next
    }
    
    # 좌표 중복 체크 + jitter
    if (any(duplicated(bayesDF[, .(lon, lat)]))) {
      cat("  [Warning]", YR, "년: 중복 좌표 감지됨 -> jitter 적용/n")
      set.seed(YR)
      bayesDF$lon <- jitter(bayesDF$lon, amount = 0.0001)
      bayesDF$lat <- jitter(bayesDF$lat, amount = 0.0001)
    }
    
    coords <- as.matrix(bayesDF[, .(lon, lat)])
    
    dist_check <- dist(coords)
    if (min(dist_check) == 0) stop("여전히 겹치는 좌표가 있습니다. Jitter 실패.")
    
    # 3. 관측값 분산 계산 (informative prior에 사용)
    S_obs2 <- var(bayesDF$month_wt_srd, na.rm = TRUE)
    
    # 4. 두 가지 prior 시나리오 정의
    # 4-1) Flat Prior: IG(0.001, 0.001)
    priors_flat <- list(
      "beta.Flat",
      "phi.Unif"    = c(0.001, 6),       # U(0.001, 6)
      "sigma.sq.IG" = c(0.001, 0.001),   # 매우 vague
      "tau.sq.IG"   = c(0.001, 0.001)
    )
    
    # 4-2) Informative Prior: IG(2, 1/S_obs^2)
    #     (설명에 나온 대로 구현)
    priors_inf <- list(
      "beta.Flat",
      "phi.Unif"    = c(0.001, 6),
      "sigma.sq.IG" = c(2, 1 / S_obs2),
      "tau.sq.IG"   = c(2, 1 / S_obs2)
    )
    
    # 5. 두 모델 피팅 (Flat / Informative)
    cat("  - MCMC 모델링 (Flat vs Informative) .../n")
    
    res_flat <- run_spLM_with_priors(
      bayesDF     = bayesDF,
      coords      = coords,
      priors      = priors_flat,
      scenario_name = "flat",
      n.samples   = n.samples,
      burn.in     = burn.in,
      starting    = starting,
      tuning      = tuning,
      cov.model   = "spherical",
      n.report    = 2000,
      seed        = YR + 100  # 시나리오별로 다른 seed 줄 수도 있음
    )
    
    res_inf <- run_spLM_with_priors(
      bayesDF     = bayesDF,
      coords      = coords,
      priors      = priors_inf,
      scenario_name = "informative",
      n.samples   = n.samples,
      burn.in     = burn.in,
      starting    = starting,
      tuning      = tuning,
      cov.model   = "spherical",
      n.report    = 2000,
      seed        = YR + 200
    )
    
    # 6. 위치별 posterior mean smoothed risk 저장
    bayesDF[, smooth_srd_flat := res_flat$beta.mean + res_flat$w.mean]
    bayesDF[, smooth_srd_inf  := res_inf$beta.mean  + res_inf$w.mean]
    
    bayesDF[, w_mean_flat := res_flat$w.mean]
    bayesDF[, w_sd_flat   := res_flat$w.sd]
    bayesDF[, w_mean_inf  := res_inf$w.mean]
    bayesDF[, w_sd_inf    := res_inf$w.sd]
    
    # 7. 연도별 결과 저장
    out_file_year <- file.path(
      OUT_BAYES_DIR,
      sprintf("%d_bayes_prior_sensitivity.csv", YR)
    )
    fwrite(bayesDF, out_file_year, bom = TRUE)
    
    cat("  >>", YR, "년 저장 완료:", out_file_year, "/n")
    cat("     · Flat  SNR (mean ± sd):",
        round(res_flat$snr.mean, 3), "±", round(res_flat$snr.sd, 3), "/n")
    cat("     · Inform SNR (mean ± sd):",
        round(res_inf$snr.mean, 3), "±", round(res_inf$snr.sd, 3), "/n")
    
    # 8. 요약 정보 저장 (연도 단위)
    summary_list[[as.character(YR)]] <- data.table(
      year                 = YR,
      S_obs2               = S_obs2,
      beta_mean_flat       = res_flat$beta.mean,
      beta_mean_inf        = res_inf$beta.mean,
      snr_mean_flat        = res_flat$snr.mean,
      snr_sd_flat          = res_flat$snr.sd,
      snr_mean_informative = res_inf$snr.mean,
      snr_sd_informative   = res_inf$snr.sd
    )
    
    # 메모리 정리
    rm(res_flat, res_inf)
    gc()
    
  }, error = function(e) {
    cat("  [Error]", YR, "년 처리 중 오류 발생:/n")
    message(e)
  })
}

cat("/n▶ 모든 연도에 대한 prior 민감도 분석 종료./n")

# ==============================================================================
# [통합 요약] 연도별 SNR 및 beta 비교 테이블 저장
# ==============================================================================

if (length(summary_list) > 0) {
  summary_dt <- rbindlist(summary_list, fill = TRUE)
  out_summary_file <- file.path(OUT_BAYES_DIR, "prior_sensitivity_summary_by_year.csv")
  fwrite(summary_dt, out_summary_file, bom = TRUE)
  cat("▶ 연도별 SNR / beta 요약 저장 완료:", out_summary_file, "/n")
} else {
  cat("※ 요약할 결과가 없습니다. (성공적으로 피팅된 연도가 없음)/n")
}

