library(data.table)
library(sp)
library(spBayes)

## ==============================================================================
## [0] 경로 및 공통 설정
## ==============================================================================

base_dir      <- "C:/Users/USER/Desktop/thesis"
OUT_OR_DIR    <- file.path(base_dir, "output3")
OUT_BAYES_DIR <- file.path(OUT_OR_DIR, "bayes_OR_prior_sensitivity")

if (!dir.exists(OUT_BAYES_DIR)) dir.create(OUT_BAYES_DIR, recursive = TRUE)

# 시군구별 OR 결과 파일 (앞에서 만든 결과)
OR_FILE   <- file.path(OUT_OR_DIR, "SGG_level_RR_2015_2024_wildfire_affected_cc_hourly.csv")

# 시군구 중심점 (좌표) 파일
cent_file <- file.path(base_dir, "SGG", "sgg_cd_2024_with_centroid.csv")

# 사용할 lag 범위
lags <- 0:1  # 0~24시간

# MCMC 설정
n.samples <- 10000
burn.in   <- 7500  # 75% burn-in (RD 민감도 분석과 동일하게)

# 초기값 및 튜닝 (두 시나리오 공통)
starting <- list("phi" = 1, "sigma.sq" = 5, "tau.sq" = 1)
tuning   <- list("phi" = 0.1, "sigma.sq" = 0.1, "tau.sq" = 0.1)


## ==============================================================================
## [1] 데이터 로드: OR + centroid
## ==============================================================================

if (!file.exists(OR_FILE)) stop("OR_FILE 이 존재하지 않습니다: ", OR_FILE)
if (!file.exists(cent_file)) stop("cent_file 이 존재하지 않습니다: ", cent_file)

or_df <- fread(OR_FILE, encoding = "UTF-8")
cent  <- fread(cent_file, encoding = "UTF-8")

# SGG_CD 문자형 통일
or_df[,  SGG_CD := as.character(SGG_CD)]
cent[,  SGG_CD := as.character(SGG_CD)]

# centroid에서 필요한 컬럼만 사용 (lon, lat 이름은 네 파일에 맞게 수정)
cent_small <- cent[, .(SGG_CD, lon, lat)]


## ==============================================================================
## [2] 헬퍼 함수: priors 바꿔가며 spLM + spRecover + SNR 계산
## ==============================================================================

run_spLM_OR_with_priors <- function(bayesDF, coords, response_col,
                                    priors, scenario_name,
                                    n.samples, burn.in,
                                    starting, tuning,
                                    cov.model = "spherical",
                                    n.report = 2000,
                                    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  cat("    · 시나리오:", scenario_name, " / n =", nrow(bayesDF), "\n")
  
  # formula 동적으로 생성: logRR ~ 1
  fmla <- as.formula(paste(response_col, "~ 1"))
  
  m.1 <- spLM(
    fmla,
    coords    = coords,
    data      = bayesDF,
    starting  = starting,
    tuning    = tuning,
    priors    = priors,
    cov.model = cov.model,
    n.samples = n.samples,
    verbose   = FALSE,
    n.report  = n.report,
    # RR 추정치의 분산을 이용한 가중치 (정밀도 기반)
    weights   = bayesDF$prec_weights
  )
  
  # 분산 파라미터 posterior
  theta.samples <- m.1$p.theta.samples
  if (burn.in >= 1) {
    theta.samples <- theta.samples[burn.in:n.samples, , drop = FALSE]
  }
  
  sigma.sq.post <- theta.samples[, "sigma.sq"]
  tau.sq.post   <- theta.samples[, "tau.sq"]
  snr.post      <- sigma.sq.post / (sigma.sq.post + tau.sq.post)  # SNR
  
  snr.mean <- mean(snr.post, na.rm = TRUE)
  snr.sd   <- sd(snr.post,  na.rm = TRUE)
  
  # 공간효과 복원
  m.1.rec <- spRecover(m.1, start = burn.in, verbose = FALSE)
  
  beta.samples <- m.1.rec$p.beta.recover.samples
  w.samples    <- m.1.rec$p.w.recover.samples
  
  beta.mean <- mean(beta.samples)
  w.mean    <- apply(w.samples, 1, mean)
  w.sd      <- apply(w.samples, 1, sd)
  
  list(
    beta.mean  = beta.mean,
    w.mean     = w.mean,
    w.sd       = w.sd,
    snr.mean   = snr.mean,
    snr.sd     = snr.sd,
    theta.post = theta.samples
  )
}


## ==============================================================================
## [3] lag별로 flat vs informative prior 민감도 분석
## ==============================================================================

summary_list <- list()  # lag별 요약 저장

for (lg in lags) {
  cat("==========================================================\n")
  cat("▶ LAG =", lg, "시간 기준 베이지안 공간 prior 민감도 분석 시작\n")
  cat("==========================================================\n")
  
  # ----------------------------------------------------------
  # 3-1. 해당 lag 데이터 필터링
  # ----------------------------------------------------------
  dat_lag <- or_df[lag_hour == lg & note == "Success"]
  dat_lag <- dat_lag[
    is.finite(RR)    & !is.na(RR)    &
      is.finite(RR_ll) & !is.na(RR_ll) &
      is.finite(RR_ul) & !is.na(RR_ul)
  ]
  dat_lag[, SGG_CD := as.character(SGG_CD)]
  
  if (nrow(dat_lag) < 10) {
    cat("  [Skip] lag", lg, ": 유효 시군구 수 부족 (n <", nrow(dat_lag), ")\n")
    next
  }
  
  # ----------------------------------------------------------
  # 3-2. logRR, SE(logRR), 분산, 정밀도(prec_weights) 계산
  # ----------------------------------------------------------
  dat_lag[, logRR    := log(RR)]
  dat_lag[, logRR_ll := log(RR_ll)]
  dat_lag[, logRR_ul := log(RR_ul)]
  
  # 95% CI 기반 SE
  dat_lag[, se_logRR := (logRR_ul - logRR_ll) / (2 * 1.96)]
  
  dat_lag[, logRR_var    := se_logRR^2]
  dat_lag[, prec_weights := 1 / logRR_var]
  
  dat_lag <- dat_lag[is.finite(prec_weights) & prec_weights > 0]
  
  if (nrow(dat_lag) < 10) {
    cat("  [Skip] lag", lg, ": 정밀도(prec) 유효 시군구 수 부족\n")
    next
  }
  
  # ----------------------------------------------------------
  # 3-3. centroid merge + 좌표 정리
  # ----------------------------------------------------------
  bayesDF <- merge(dat_lag, cent_small, by = "SGG_CD", all.x = TRUE)
  bayesDF <- bayesDF[!is.na(lon) & !is.na(lat)]
  
  if (nrow(bayesDF) < 10) {
    cat("  [Skip] lag", lg, ": 좌표가 있는 시군구 n <", nrow(bayesDF), "\n")
    next
  }
  
  # 좌표 중복 시 jitter
  if (any(duplicated(bayesDF[, .(lon, lat)]))) {
    cat("  [Warning] lag", lg, ": 중복 좌표 감지됨 -> jitter 적용\n")
    set.seed(1000 + lg)
    bayesDF$lon <- jitter(bayesDF$lon, amount = 0.0001)
    bayesDF$lat <- jitter(bayesDF$lat, amount = 0.0001)
  }
  coords <- as.matrix(bayesDF[, .(lon, lat)])
  
  d_check <- dist(coords)
  if (min(d_check) == 0) {
    cat("  [Error] lag", lg, ": jitter 후에도 겹치는 좌표 존재 → 스킵\n")
    next
  }
  
  # ----------------------------------------------------------
  # 3-4. 관측 logRR 분산 계산 (informative prior에 사용)
  # ----------------------------------------------------------
  S_obs2 <- var(bayesDF$logRR, na.rm = TRUE)
  if (!is.finite(S_obs2) || S_obs2 <= 0) {
    cat("  [Skip] lag", lg, ": 관측 logRR 분산 S_obs2 비정상\n")
    next
  }
  
  # ----------------------------------------------------------
  # 3-5. 두 종류 prior 정의
  #   - priors_flat : 매우 vague IG (flat에 가까운)
  #   - priors_inf  : IG(2, 1/S_obs2) 형태의 정보 있는 prior
  # ----------------------------------------------------------
  
  # (1) Flat-ish prior (매우 약한 prior)
  priors_flat <- list(
    "beta.Flat",
    "phi.Unif"    = c(0.001, 6),      # 공간 범위: 거의 non-informative
    "sigma.sq.IG" = c(0.001, 0.001),  # 아주 퍼진 IG
    "tau.sq.IG"   = c(0.001, 0.001)
  )
  
  # (2) Informative prior: 데이터 분산 S_obs2를 이용
  priors_inf <- list(
    "beta.Flat",
    "phi.Unif"    = c(0.001, 6),
    "sigma.sq.IG" = c(2, 1 / S_obs2),
    "tau.sq.IG"   = c(2, 1 / S_obs2)
  )
  
  # ----------------------------------------------------------
  # 3-6. spLM + spRecover (Flat vs Informative)
  # ----------------------------------------------------------
  cat("  - MCMC (Flat vs Informative) 실행 중...\n")
  
  res_flat <- run_spLM_OR_with_priors(
    bayesDF       = bayesDF,
    coords        = coords,
    response_col  = "logRR",
    priors        = priors_flat,
    scenario_name = "flat",
    n.samples     = n.samples,
    burn.in       = burn.in,
    starting      = starting,
    tuning        = tuning,
    cov.model     = "spherical",
    n.report      = 2000,
    seed          = 1000 + lg
  )
  
  res_inf <- run_spLM_OR_with_priors(
    bayesDF       = bayesDF,
    coords        = coords,
    response_col  = "logRR",
    priors        = priors_inf,
    scenario_name = "informative",
    n.samples     = n.samples,
    burn.in       = burn.in,
    starting      = starting,
    tuning        = tuning,
    cov.model     = "spherical",
    n.report      = 2000,
    seed          = 2000 + lg
  )
  
  # ----------------------------------------------------------
  # 3-7. 위치별 smooth logRR / RR 계산 및 저장
  # ----------------------------------------------------------
  bayesDF[, smooth_logRR_flat := res_flat$beta.mean + res_flat$w.mean]
  bayesDF[, smooth_RR_flat    := exp(smooth_logRR_flat)]
  bayesDF[, smooth_logRR_inf  := res_inf$beta.mean  + res_inf$w.mean]
  bayesDF[, smooth_RR_inf     := exp(smooth_logRR_inf)]
  
  bayesDF[, w_mean_flat := res_flat$w.mean]
  bayesDF[, w_sd_flat   := res_flat$w.sd]
  bayesDF[, w_mean_inf  := res_inf$w.mean]
  bayesDF[, w_sd_inf    := res_inf$w.sd]
  
  bayesDF[, lag_hour := lg]
  
  out_file_lag <- file.path(
    OUT_BAYES_DIR,
    sprintf("lag%02d_bayes_OR_prior_sensitivity.csv", lg)
  )
  fwrite(bayesDF, out_file_lag, bom = TRUE)
  
  cat("  >> lag", lg, "저장 완료:", out_file_lag, "\n")
  cat("     · Flat  SNR (mean ± sd):",
      round(res_flat$snr.mean, 3), "±", round(res_flat$snr.sd, 3), "\n")
  cat("     · Inform SNR (mean ± sd):",
      round(res_inf$snr.mean, 3), "±", round(res_inf$snr.sd, 3), "\n")
  
  # ----------------------------------------------------------
  # 3-8. lag 단위 요약 저장
  # ----------------------------------------------------------
  summary_list[[as.character(lg)]] <- data.table(
    lag_hour             = lg,
    n_region             = nrow(bayesDF),
    S_obs2               = S_obs2,
    beta_mean_flat       = res_flat$beta.mean,
    beta_mean_inf        = res_inf$beta.mean,
    snr_mean_flat        = res_flat$snr.mean,
    snr_sd_flat          = res_flat$snr.sd,
    snr_mean_informative = res_inf$snr.mean,
    snr_sd_informative   = res_inf$snr.sd
  )
  
  rm(res_flat, res_inf)
  gc()
}

cat("\n▶ 모든 lag에 대한 prior 민감도 분석 종료.\n")

## ==============================================================================
## [4] lag별 SNR / beta 요약 테이블 저장
## ==============================================================================

if (length(summary_list) > 0) {
  summary_dt <- rbindlist(summary_list, fill = TRUE)
  out_summary_file <- file.path(OUT_BAYES_DIR, "prior_sensitivity_summary_by_lag.csv")
  fwrite(summary_dt, out_summary_file, bom = TRUE)
  cat("▶ lag별 SNR / beta 요약 저장 완료:", out_summary_file, "\n")
} else {
  cat("※ 요약할 결과가 없습니다. (성공적으로 피팅된 lag 없음)\n")
}
