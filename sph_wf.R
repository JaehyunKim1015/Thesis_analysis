## ================================================================
## 0. 패키지 로드
## ================================================================
library(data.table)
library(sp)
library(spBayes)

## ================================================================
## 1. 경로 / 파일 설정
## ================================================================
base_dir    <- "C:/Users/USER/Desktop/thesis"

# hourly CC 결과 (시군구별 RR, RR_ll, RR_ul, lag_hour, note 등)
OR_FILE     <- file.path(base_dir, "output3",
                         "SGG_level_RR_2015_2024_wildfire_affected_cc_hourly.csv")
# 시군구 centroid (lon, lat 포함)
CENT_FILE   <- file.path(base_dir, "SGG", "sgg_cd_2024_with_centroid.csv")

# 베이지안 공간모형 결과 저장 폴더
OUT_OR_DIR  <- file.path(base_dir, "output3", "bayes_OR")
if (!dir.exists(OUT_OR_DIR)) dir.create(OUT_OR_DIR, recursive = TRUE)

# 분석할 lag (시간 단위)
lags <- c(0, 1)

## ================================================================
## 2. 데이터 로드
## ================================================================
if (!file.exists(OR_FILE)) stop("OR_FILE 이 존재하지 않습니다: ", OR_FILE)
if (!file.exists(CENT_FILE)) stop("CENT_FILE 이 존재하지 않습니다: ", CENT_FILE)

or_df <- fread(OR_FILE, encoding = "UTF-8")
cent  <- fread(CENT_FILE, encoding = "UTF-8")

# SGG_CD를 문자형으로 맞추기
or_df[,  SGG_CD := as.character(SGG_CD)]
cent[,  SGG_CD := as.character(SGG_CD)]

# centroid 파일에서 사용할 컬럼만 남기기 (lon/lat 컬럼 이름은 파일에 맞게 확인)
# 예: SGG_CD, lon, lat 가 있다고 가정
cent <- cent[, .(SGG_CD, lon, lat)]

## ================================================================
## 3. spLM 설정 (Gaussian 공간 모형)
##    - logRR ~ 1 + spatial random effect
## ================================================================
n.samples <- 10000
burn.in   <- 0.5 * n.samples   # 50% burn-in

starting <- list(
  "phi"      = 1,   # 공간 decay 초기값
  "sigma.sq" = 0.5, # 공간 분산
  "tau.sq"   = 0.5  # nugget (비공간 잡음분산)
)

tuning <- list(
  "phi"      = 0.1,
  "sigma.sq" = 0.1,
  "tau.sq"   = 0.1
)

priors <- list(
  "beta.Flat",
  "phi.Unif"    = c(0.01, 10),
  "sigma.sq.IG" = c(2, 1),
  "tau.sq.IG"   = c(2, 1)
)

## ================================================================
## 4. lag별로 logRR / se_logRR 계산 + spLM 베이지안 공간모형
## ================================================================
all_bayes_results <- list()

for (lg in lags) {
  
  cat("==========================================================\n")
  cat("▶ LAG =", lg, "시간 기준 베이지안 공간 분석 시작\n")
  cat("==========================================================\n")
  
  ## ---------------------------
  ## 4-1. 해당 lag 데이터 필터링
  ## ---------------------------
  dat_lag <- or_df[lag_hour == lg & note == "Success"]
  dat_lag <- dat_lag[
    is.finite(RR)    & !is.na(RR) &
      is.finite(RR_ll) & !is.na(RR_ll) &
      is.finite(RR_ul) & !is.na(RR_ul)
  ]
  dat_lag[, SGG_CD := as.character(SGG_CD)]
  
  if (nrow(dat_lag) < 10) {
    cat("  [Skip] lag", lg, "유효 시군구 수 부족 (n <", 10, ")\n\n")
    next
  }
  
  ## ---------------------------
  ## 4-2. logRR / SE / Var 계산
  ## ---------------------------
  dat_lag[, logRR    := log(RR)]
  dat_lag[, logRR_ll := log(RR_ll)]
  dat_lag[, logRR_ul := log(RR_ul)]
  
  # 95% CI 를 이용한 SE(logRR) 추정
  dat_lag[, se_logRR := (logRR_ul - logRR_ll) / (2 * 1.96)]
  dat_lag[, logRR_var := se_logRR^2]
  
  # 정밀도(precision) = 1 / Var (meta-analysis에서 자주 쓰는 개념)
  dat_lag[, prec_weights := 1 / logRR_var]
  
  # 아주 이상한 값(무한대, 0, NA 등) 제거
  dat_lag <- dat_lag[is.finite(prec_weights) & prec_weights > 0]
  
  if (nrow(dat_lag) < 10) {
    cat("  [Skip] lag", lg, "에서 정밀도(prec)가 유효한 시군구 수 부족.\n\n")
    next
  }
  
  ## ---------------------------
  ## 4-3. centroid 병합 + 좌표 jitter
  ## ---------------------------
  bayesDF <- merge(dat_lag, cent, by = "SGG_CD", all.x = TRUE)
  bayesDF <- bayesDF[!is.na(lon) & !is.na(lat)]
  
  if (nrow(bayesDF) < 10) {
    cat("  [Skip] lag", lg, " centroid 병합 후 유효 시군구 수 부족.\n\n")
    next
  }
  
  # 같은 좌표(겹치는 시군구) 있으면 조금 흔들어서(dpotrf 에러 방지)
  if (any(duplicated(bayesDF[, .(lon, lat)]))) {
    cat("  - 중복 좌표 감지 → jitter 적용\n")
    set.seed(1000 + lg)
    bayesDF$lon <- jitter(bayesDF$lon, amount = 0.0001)
    bayesDF$lat <- jitter(bayesDF$lat, amount = 0.0001)
  }
  coords <- as.matrix(bayesDF[, .(lon, lat)])
  
  # 거리 0 인 지점이 남아있는지 마지막 체크
  if (min(dist(coords)) == 0) {
    stop("여전히 동일한 좌표가 있습니다. jitter amount를 더 크게 해 주세요.")
  }
  
  ## ---------------------------
  ## 4-4. spLM 실행 (Gaussian)
  ##    logRR ~ 1 + 공간적 random effect
  ##    ⚠ spLM 은 observation-specific weights 인자를 직접 받지 못함.
  ##      (prec_weights는 나중에 참고용으로만 사용 / INLA 등에서 활용 가능)
  ## ---------------------------
  cat("  - spLM MCMC 실행 중... (시군구 수:", nrow(bayesDF), ")\n")
  
  m.1 <- spLM(
    logRR ~ 1,
    coords    = coords,
    data      = bayesDF,
    starting  = starting,
    tuning    = tuning,
    priors    = priors,
    cov.model = "spherical",
    n.samples = n.samples,
    verbose   = FALSE,
    n.report  = 2000
  )
  
  ## ---------------------------
  ## 4-5. Posterior 복원
  ## ---------------------------
  cat("  - spRecover (burn-in 이후 샘플 사용) ...\n")
  m.1.rec <- spRecover(m.1, start = burn.in, verbose = FALSE)
  
  beta.samples <- m.1.rec$p.beta.recover.samples    # 고정효과(전국 평균 logRR)
  w.samples    <- m.1.rec$p.w.recover.samples       # 공간 random effect
  
  beta_mean <- mean(beta.samples)
  w_mean    <- apply(w.samples, 1, mean)
  w_sd      <- apply(w.samples, 1, sd)
  
  ## ---------------------------
  ## 4-6. smooth logRR / RR 계산
  ## ---------------------------
  bayesDF[, beta_mean    := beta_mean]
  bayesDF[, w_mean       := w_mean]
  bayesDF[, w_sd         := w_sd]
  bayesDF[, smooth_logRR := beta_mean + w_mean]
  bayesDF[, smooth_RR    := exp(smooth_logRR)]
  
  ## ---------------------------
  ## 4-7. 저장 및 통합
  ## ---------------------------
  out_file <- file.path(OUT_OR_DIR,
                        sprintf("lag%02d_bayes_OR_result.csv", lg))
  fwrite(bayesDF, out_file, bom = TRUE)
  
  cat("  >> lag", lg, "완료! 평균 Log(RR):", round(beta_mean, 4), "\n")
  cat("  >> 결과 저장:", out_file, "\n\n")
  
  # 후속 통합 분석용으로 필요한 것만 저장
  bayesDF[, lag_hour := lg]
  all_bayes_results[[as.character(lg)]] <-
    bayesDF[, .(SGG_CD,
                lag_hour,
                smooth_RR,
                smooth_logRR,
                logRR,
                logRR_var,
                prec_weights,
                n_case)]
  
  # 메모리 정리
  rm(m.1, m.1.rec, beta.samples, w.samples)
  gc()
}

cat("\n▶ 모든 lag에 대한 베이지안 공간 분석 종료.\n")

## ================================================================
## 5. lag 통합 결과 한 번에 저장 (선택사항)
## ================================================================
if (length(all_bayes_results) > 0) {
  all_bayes_df <- rbindlist(all_bayes_results, use.names = TRUE, fill = TRUE)
  out_all_file <- file.path(OUT_OR_DIR, "all_lags_bayes_OR_result_combined.csv")
  fwrite(all_bayes_df, out_all_file, bom = TRUE)
  cat("▶ 모든 lag 결과 통합 파일 저장 완료:", out_all_file, "\n")
} else {
  cat("⚠ 유효한 베이지안 결과가 없습니다. (모든 lag에서 실패했을 수 있음)\n")
}
