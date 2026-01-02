library(data.table)
library(lubridate)
library(survival)
library(splines)
library(dlnm)    

## =========================================
## 0. 경로 및 연도 설정
## =========================================
base_dir   <- "C:/Users/USER/Desktop/thesis"
expo_dir   <- file.path(base_dir, "wf_bi")
nedis_dir  <- file.path(base_dir, "NEDIS_count_daily")
wt_dir     <- file.path(base_dir,"wildfire","wt_daily")
# SGG_CD 매핑 경로 (sgg_dir) 제거됨
out_dir    <- file.path(base_dir, "output1")
dir.create(out_dir, showWarnings = FALSE)

years <- 2015:2024

expo_list <- list()
ha_list   <- list()
wt_list   <- list()

## =========================================
## 1. 연도별 노출 / 건강 / 날씨 데이터 읽기
## =========================================
for (year_use in years) {
  cat("\n===============================\n")
  cat(" Reading year:", year_use, "\n")
  cat("===============================\n")
  
  ## -----------------------------
  ## 1-1. Exposure (연도별_wfbi.csv)
  ## -----------------------------
  expo_file <- file.path(expo_dir, paste0(year_use, "_wfbi_with_sgg.csv"))
  cat("  expo_file:", expo_file, "\n")
  expo <- fread(expo_file)
  
  cat("  raw nrow(expo):", nrow(expo), "\n")
  
  ## SGG 코드 이름 통일 및 타입 변환
  sgg_col_e <- grep("SGG_CD", names(expo), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(sgg_col_e)) stop("SGG column not found in expo for year ", year_use)
  setnames(expo, old = sgg_col_e, new = "SGG_CD")
  expo[, SGG_CD := as.character(SGG_CD)]
  
  ## 날짜 형식 변환 및 컬럼 이름 통일
  expo[, date := ymd(date)]
  
  ## 해당 연도만 필터
  expo <- expo[year(date) == year_use]
  
  ## wildfire binary → wf 로 통일 (smoke_day 사용)
  wf_col <- grep("smoke_day|wf", names(expo), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(wf_col)) stop("wf column (smoke_day) not found in expo for year ", year_use)
  setnames(expo, old = wf_col, new = "wf")
  
  expo[, year := year_use]
  expo_list[[as.character(year_use)]] <- expo
  
  ## -----------------------------
  ## 1-2. Health (NEDIS_count)
  ## -----------------------------
  ha_file <- file.path(nedis_dir, paste0(year_use, "_daily_with_centroid.csv"))
  cat("  ha_file  :", ha_file, "\n")
  ha <- fread(ha_file)
  cat("  raw nrow(ha):", nrow(ha), "\n")
  
  ## SGG 코드 통일
  sgg_col_h <- grep("SGG_CD", names(ha), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(sgg_col_h)) stop("SGG column not found in ha for year ", year_use)
  setnames(ha, old = sgg_col_h, new = "SGG_CD")
  ha[, SGG_CD := as.character(SGG_CD)]
  
  ## 날짜 컬럼 감지 및 'case_date'로 이름 통일
  date_col_h <- grep("^date$", names(ha), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(date_col_h)) stop("Date column not found in ha for year ", year_use)
  setnames(ha, old = date_col_h, new = "case_date")
  
  ## count 컬럼 이름을 'respiratory'로 통일
  count_col_h <- grep("daily_count|respiratory", names(ha), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(count_col_h) && count_col_h != "respiratory") {
    setnames(ha, old = count_col_h, new = "respiratory")
  }
  
  ## 날짜 형식 변환 (ymd()로 오류 방지)
  ha[, case_date := ymd(case_date)]
  
  if (!is.Date(ha$case_date)) stop("Date conversion failed in ha for year ", year_use)
  
  ## 해당 연도 필터
  ha <- ha[year(case_date) == year_use]
  ha[, year := year_use]
  
  ## time-stratified 변수
  ha[, `:=`(
    wday  = wday(case_date),
    month = month(case_date)
  )]
  
  ha_list[[as.character(year_use)]] <- ha
  
  ## -----------------------------
  ## 1-3. Weather (wt_daily, temp만 사용)
  ## -----------------------------
  wt_file <- file.path(wt_dir, paste0(year_use, "_wtd.csv"))
  cat("  wt_file  :", wt_file, "\n")
  wt <- fread(wt_file)
  
  wt[, date := ymd(date)]
  # temp만 사용
  wt <- wt[, .(date, SGG_CD, temp)]
  wt[, SGG_CD := as.character(SGG_CD)]
  wt[, year := year_use]
  
  wt_list[[as.character(year_use)]] <- wt
}

## 합치기
expo_all <- rbindlist(expo_list, use.names = TRUE, fill = TRUE)
ha_all   <- rbindlist(ha_list,   use.names = TRUE, fill = TRUE)
wt_all   <- rbindlist(wt_list,   use.names = TRUE, fill = TRUE)

cat("\nTotal exposure rows:", nrow(expo_all), "\n")
cat("Total health rows   :", nrow(ha_all), "\n")
cat("Total weather rows  :", nrow(wt_all), "\n")


## =========================================
## 2. time-stratified 변수 만들기
## (기존 Section 2 제거 후, Section 3이 2로 변경)
## =========================================
ha_all[, `:=`(
  wday  = wday(case_date),
  month = month(case_date),
  year  = year(case_date)
)]

expo_all[, `:=`(
  wday  = wday(date),
  month = month(date),
  year  = year(date)
)]

wt_all[, `:=`(
  wday  = wday(date),
  month = month(date),
  year  = year(date)
)]


## =========================================
## 3. case-crossover용 dt 만들기
## =========================================

## 사건날만 고르기
resp <- ha_all[respiratory > 0]
resp[, id := .I] # 각 사건(case)에 고유 ID 부여

## expo + weather 먼저 merge (SGG_CD + date 기준)
expo_wt <- merge(
  expo_all,
  wt_all,
  by = c("SGG_CD", "date", "year", "wday", "month"),
  all.x = TRUE
)

## case-crossover merge (time-stratified: SGG_CD + wday + month + year)
dt <- merge(
  resp,
  expo_wt,
  by = c("SGG_CD", "wday", "month", "year"),
  all.x = TRUE,
  allow.cartesian = TRUE
)

## case indicator (매칭된 행이 사건날이면 1, 아니면 0)
dt[, case := fifelse(case_date == date, 1L, 0L)]
setnames(dt, "respiratory", "outcome")


## =========================================
## 4. SGG 코드 유의성 및 기여도 확인 (Case Contribution Check)
## (기존 Section 4가 4로 유지)
## =========================================
sgg_contribution <- dt[case == 1, .N, by = SGG_CD][order(-N)]

cat("\n=== 분석에 포함된 SGG 코드 기여도 (Case 발생 기준) ===\n")
cat("총 유의미하게 기여하는 SGG_CD 개수 (case > 0):", nrow(sgg_contribution), "\n")
cat("가장 많이 기여한 상위 10개 SGG_CD:\n")
print(sgg_contribution[1:10])

fwrite(sgg_contribution, file.path(out_dir, "sgg_case_contribution.csv"))
cat(sprintf("-> SGG별 기여도 목록이 %s에 저장되었습니다.\n", 
            file.path(out_dir, "sgg_case_contribution.csv")))


## =========================================
## 5. DLNM용 노출 변수 설정 및 필터링
## (기존 Section 5가 5로 유지)
## =========================================
if (!"PM25_wildfire_excess" %in% names(dt)) {
  stop("PM25_wildfire_excess not found in dt")
}
dt[, pm := PM25_wildfire_excess]

## temp NA 처리 
dt <- dt[!is.na(pm) & !is.na(temp)]

cat("nrow(dt) after filtering:", nrow(dt), "\n")
cat("summary(pm):\n"); print(summary(dt$pm))
cat("temp NA proportion:", mean(is.na(dt$temp)), "\n")


## =========================================
## 6. DLNM cross-basis (lag 0–6)
## (기존 Section 6이 6으로 변경)
## =========================================
cb_pm <- crossbasis(
  dt$pm,
  lag   = 6,               # 0~6일
  argvar = list(fun = "lin"),    # PM 선형
  arglag = list(fun = "ns", df = 3)  # lag 방향 spline
)


## =========================================
## 7. clogit 모델 및 RR 계산
## (기존 Section 7/8이 7로 통합)
## =========================================

## (A) temp linear
m_dlnm_temp <- clogit(
  case ~ cb_pm + temp + strata(id),
  data    = dt,
  weights = outcome,
  method  = "approximate"
)

cat("\n=== DLNM + temp(linear) 모델 요약 ===\n")
print(summary(m_dlnm_temp))

## 누적 RR (cumulative RR, lag 0–6)
cp_cum <- crosspred(
  cb_pm,
  model = m_dlnm_temp,
  cen   = 0,
  at    = 10,
  cumul = TRUE
)

cRR <- cp_cum$cumfit

cat("\n=== 누적 RR (0→10 μg/m³, lag0–6) ===\n")
print(cRR)

## lag별 RR (lag 0–6) 추출 및 저장
cp_lag <- crosspred(cb_pm, model=m_dlnm_temp, cen=0, at=10, cumul=FALSE)
logRR    <- as.numeric(cp_lag$matfit)
lags_vec <- 0:(length(logRR) - 1)

lag_rr_dt <- data.table(
  lag      = lags_vec,
  logRR    = logRR,
  RR       = exp(logRR)
)

cat("\n=== lag별 RR (0–6, at = 10 μg/m³ 증가) ===\n")
print(lag_rr_dt)

# 결과 파일 저장 (덮어쓰지 않음)
fwrite(lag_rr_dt,
       file.path(out_dir, "dlnm_temp_only_lag0_6_rr_final.csv"))
cat(sprintf("-> lag별 RR 결과가 %s에 저장되었습니다.\n", 
            file.path(out_dir, "dlnm_temp_only_lag0_6_rr_final.csv")))