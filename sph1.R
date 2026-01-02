library(data.table)
library(sp)
library(spBayes)
library(coda)

# ==============================================================================
# [설정] 경로 및 파라미터
# ==============================================================================

OUT_RD_DIR    <- "C:/Users/USER/Desktop/thesis/output1/sgg_rd"
OUT_BAYES_DIR <- file.path(OUT_RD_DIR, "bayes_result")

if (!dir.exists(OUT_BAYES_DIR)) dir.create(OUT_BAYES_DIR, recursive = TRUE)

years <- 2015:2024

# MCMC 설정
n.samples <- 10000
burn.in   <- 0.75 * n.samples 

# 초기값 및 튜닝
starting <- list("phi" = 1, "sigma.sq" = 5, "tau.sq" = 1)
tuning   <- list("phi" = 0.1, "sigma.sq" = 0.1, "tau.sq" = 0.1)

# 사전분포 (Priors)
priors <- list(
  "beta.Flat",
  "phi.Unif" = c(0.01, 10),      # 공간 범위 (하한선을 0.01로 낮춤)
  "sigma.sq.IG" = c(2, 5),      
  "tau.sq.IG" = c(2, 1)         
)

# ==============================================================================
# [루프 실행] 2015 ~ 2023
# ==============================================================================

for (YR in years) {
  cat("==========================================================/n")
  cat("▶", YR, "년 베이지안 공간 분석 시작/n")
  cat("==========================================================/n")
  
  tryCatch({
    
    # 1. 데이터 로드
    file_path <- file.path(OUT_RD_DIR, sprintf("%d_sgg_rd_monthwt.csv", YR))
    
    if (!file.exists(file_path)) {
      cat("  [Skip]", YR, "년 입력 파일이 없습니다./n")
      next 
    }
    
    # 데이터 읽기
    baz <- fread(file_path, encoding = "UTF-8") # 혹은 "Unknown"
    
    # 2. 전처리
    bayesDF <- baz[, .(SGG_CD, month_wt_srd, lon, lat, pop)]
    bayesDF <- bayesDF[!is.na(month_wt_srd) & is.finite(month_wt_srd) & 
                         !is.na(lon) & !is.na(lat) & !is.na(pop)]
    
    # 인구 필터링
    bayesDF <- bayesDF[pop > 1000]
    
    if (nrow(bayesDF) < 10) {
      cat("  [Skip] 유효 데이터 부족./n")
      next
    }
    
    # ★★★ [핵심 수정] 중복 좌표 해결 (Jittering) ★★★
    # 좌표가 완전히 같은 지점이 하나라도 있으면 dpotrf 에러 발생함.
    # 아주 미세한 값(amount=0.0001, 약 10m)을 더해서 겹침 방지
    if (any(duplicated(bayesDF[, .(lon, lat)]))) {
      cat("  [Warning] 중복 좌표 감지됨 -> 미세 조정(Jitter) 적용/n")
      set.seed(YR) 
      bayesDF$lon <- jitter(bayesDF$lon, amount = 0.0001)
      bayesDF$lat <- jitter(bayesDF$lat, amount = 0.0001)
    }
    
    # 3. spLM 모델링
    coords <- as.matrix(bayesDF[, .(lon, lat)])
    
    cat("  - MCMC 모델링 중... (", nrow(bayesDF), "개 지역 )/n")
    
    # 에러가 자주 나는 구간이므로 안전장치 한번 더 확인
    check_dist <- dist(coords)
    if(min(check_dist) == 0) stop("여전히 겹치는 좌표가 있습니다. Jitter가 작동하지 않았습니다.")
    
    m.1 <- spLM(month_wt_srd ~ 1, 
                coords = coords,
                data = bayesDF, 
                starting = starting,
                tuning = tuning,
                priors = priors,
                cov.model = "spherical",
                n.samples = n.samples,
                verbose = FALSE,
                n.report = 2000)
    
    # 4. 공간 효과 복원 (Recover)
    cat("  - 공간 효과 복원 중.../n")
    m.1.rec <- spRecover(m.1, start = burn.in, verbose = FALSE)
    
    # 결과 추출
    beta.samples <- m.1.rec$p.beta.recover.samples
    w.samples    <- m.1.rec$p.w.recover.samples
    
    # 5. 통계치 계산
    beta_mean <- mean(beta.samples) 
    w_mean    <- apply(w.samples, 1, mean) 
    w_sd      <- apply(w.samples, 1, sd)
    
    # 최종 보정값
    bayesDF$smooth_srd <- beta_mean + w_mean
    bayesDF$w_mean <- w_mean
    bayesDF$w_sd   <- w_sd
    
    # 6. 저장
    out_file <- file.path(OUT_BAYES_DIR, sprintf("%d_bayes_result.csv", YR))
    fwrite(bayesDF, out_file, bom = TRUE)
    
    cat("  >> 완료! 전국 평균(Beta):", round(beta_mean, 4), "/n")
    
    # 메모리 정리
    rm(m.1, m.1.rec, beta.samples, w.samples)
    gc() 
    
  }, error = function(e) {
    cat("  [Error]", YR, "년 처리 중 오류 발생:/n")
    message(e)
  })
}

cat("/n▶ 모든 연도 작업 종료./n")

library(data.table)

# ==============================================================================
# 7. 9년치 결과 통합 및 행정구역 변경 검증
# ==============================================================================

# 베이지안 결과가 저장된 폴더
OUT_BAYES_DIR <- "C:/Users/USER/Desktop/thesis/output1/sgg_rd/bayes_result"

# 파일 목록 가져오기
result_files <- list.files(OUT_BAYES_DIR, pattern = "_bayes_result.csv", full.names = TRUE)

if (length(result_files) == 0) {
  stop("저장된 베이지안 결과 파일이 없습니다.")
}

# 1. 파일 하나로 합치기
all_results <- lapply(result_files, function(f) {
  dt <- fread(f, encoding = "UTF-8")
  # 파일명에서 연도 추출
  yr <- as.numeric(gsub(".*([0-9]{4})_bayes_result.csv", "//1", f))
  dt[, year := yr]
  
  # 필요한 컬럼만 선택 (SGG_CD가 핵심)
  return(dt[, .(SGG_CD, year, smooth_srd, w_mean, w_sd, lon, lat, pop)])
})

final_df <- rbindlist(all_results, fill = TRUE)

# -------------------------------------------------------
# 2. [검증] 시군구 코드 생존 기간 확인 (문제 지역 찾기)
# -------------------------------------------------------
sgg_check <- final_df[, .(
  years_count = .N,               # 데이터가 존재하는 연도 수 (최대 9)
  min_year    = min(year),        # 시작 연도
  max_year    = max(year),        # 끝 연도
  mean_risk   = mean(smooth_srd, na.rm = TRUE)
), by = SGG_CD]

# 문제가 있는(9년치가 다 없는) 지역들 추출
problem_sgg <- sgg_check[years_count < 9]

cat("======================================================/n")
cat("▶ 행정구역 변경 의심 지역 분석/n")
cat("======================================================/n")
cat("  - 전체 시군구 코드 개수:", nrow(sgg_check), "/n")
cat("  - 9년 개근(정상) 시군구:", nrow(sgg_check[years_count == 9]), "/n")
cat("  - 중간에 변경된 시군구:", nrow(problem_sgg), "/n")

if (nrow(problem_sgg) > 0) {
  cat("/n[주의] 아래 시군구들은 9년치 데이터가 온전히 없습니다./n")
  print(problem_sgg[order(SGG_CD)])
  
  # 별도 파일로 저장해서 눈으로 확인해보세요
  fwrite(problem_sgg, file.path(OUT_BAYES_DIR, "sgg_code_changes_check.csv"), bom = TRUE)
  cat("  >> 변경 이력 상세 파일 저장됨: sgg_code_changes_check.csv/n")
} else {
  cat("/n[완벽] 모든 시군구가 9년 내내 유지되었습니다! (매우 드문 경우입니다)/n")
}

# -------------------------------------------------------
# 3. 최종 평균 계산 (9년치 통합)
# -------------------------------------------------------
# 변경된 지역이라도 '존재했던 기간 동안의 평균'을 구해서 지도에 표출하는 것이 일반적입니다.
sgg_9year_avg <- final_df[, .(
  smooth_srd = mean(smooth_srd, na.rm = TRUE),
  w_mean     = mean(w_mean, na.rm = TRUE),
  lon        = mean(lon, na.rm = TRUE), # 좌표가 미세하게 바뀌었을 수 있으므로 평균 사용
  lat        = mean(lat, na.rm = TRUE),
  data_years = .N  # 몇 년치 데이터 기반인지 표시 (나중에 지도 그릴 때 참고)
), by = SGG_CD]

# 저장
fwrite(final_df, file.path(OUT_BAYES_DIR, "all_years_bayes_combined.csv"), bom = TRUE)
fwrite(sgg_9year_avg, file.path(OUT_BAYES_DIR, "sgg_9year_average.csv"), bom = TRUE)

cat("/n▶ 통합 저장 완료./n")