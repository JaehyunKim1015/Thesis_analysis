library(data.table)
library(lubridate)
library(survival)
library(splines)
library(dlnm)

## =========================================
## 0. 경로 및 연도 설정
## =========================================
base_dir  <- "C:/Users/USER/Desktop/thesis"

expo_dir  <- file.path(base_dir, "wildfire", "total_clean")      # ✅ hourly total_clean
nedis_dir <- file.path(base_dir, "NEDIS_count_hourly")           # ✅ hourly NEDIS
map_dir   <- file.path(base_dir, "sgg_cd")
out_dir   <- file.path(base_dir, "output3")                      # 결과 저장 폴더 (hourly용)
dir.create(out_dir, showWarnings = FALSE)

years <- 2015:2024

## ==========================================
## 0.5. SGG 매핑 테이블 생성 (2024 기준)
## ==========================================
cat("Generating robust SGG mapping table (Target 2024 standard, all exceptions)...\n")

map_file_2022 <- file.path(map_dir, "sgg_cd_2022.csv")
map_file_2023 <- file.path(map_dir, "sgg_cd_2023.csv")
map_file_2024 <- file.path(map_dir, "sgg_cd_2024.csv")

if (!all(file.exists(map_file_2022, map_file_2023, map_file_2024))) {
  stop(paste("❌ Error: SGG mapping files not found in", map_dir))
}

# 1. 데이터 로드 및 전처리
df_2022 <- fread(map_file_2022)
df_2023 <- fread(map_file_2023)
df_2024 <- fread(map_file_2024)

setnames(df_2022, c("code", "sido", "sigungu"))
setnames(df_2023, c("code", "sido", "sigungu"))
setnames(df_2024, c("code", "sido", "sigungu"))

df_2022[, code := as.character(code)]
df_2023[, code := as.character(code)]
df_2024[, code := as.character(code)]

remove_space <- function(x) gsub("\\s+", "", x)
df_2022[, `:=`(sido = remove_space(sido), sigungu = remove_space(sigungu))]
df_2023[, `:=`(sido = remove_space(sido), sigungu = remove_space(sigungu))]
df_2024[, `:=`(sido = remove_space(sido), sigungu = remove_space(sigungu))]

# 2. Historical 코드 통합
df_hist <- rbindlist(list(
  df_2022[, .(old_code = code, sido, sigungu)],
  df_2023[, .(old_code = code, sido, sigungu)]
), use.names = TRUE)
df_hist <- unique(df_hist)

# 3. Target List (2024 기준)
target_list <- unique(df_2024[, .(SGG_CD = code)])

# 4. Base Mapping (sido + sigungu 기준)
base_map <- merge(
  df_hist,
  df_2024[, .(new_code = code, sido, sigungu)],
  by = c("sido", "sigungu"),
  all.x = TRUE
)

# 5. 예외 처리
# 5-1. 군위군: 경북(47720) -> 대구(27720)
base_map[sigungu == "군위군" & old_code == "47720" & is.na(new_code), new_code := "27720"] 

# 5-2. 인천 남구 -> 미추홀구
michuhol_code_2024 <- df_2024[sigungu %like% "미추홀", code][1]
if (!is.na(michuhol_code_2024)) {
  base_map[old_code == "28170" & is.na(new_code), new_code := michuhol_code_2024]
}

# 5-3. 부천시 (3구 → 통합 → 3구)
bucheon_target_2024 <- df_2024[sido == "경기" & sigungu %like% "원미구", code][1]

if (!is.na(bucheon_target_2024)) {
  old_3gu_codes <- c("41195", "41197", "41199")
  base_map[old_code %in% old_3gu_codes & is.na(new_code), new_code := bucheon_target_2024]
  
  bucheon_single_codes <- df_hist[sigungu == "부천시" & sido == "경기", old_code]
  base_map[old_code %in% bucheon_single_codes & is.na(new_code), new_code := bucheon_target_2024]
}

# 6. 최종 매핑 테이블
final_map_list <- list()
final_map_list[[1]] <- base_map[!is.na(old_code) & !is.na(new_code), .(old_code, new_code)]
existing_olds <- final_map_list[[1]]$old_code
self_map <- target_list[!SGG_CD %in% existing_olds, .(old_code = SGG_CD, new_code = SGG_CD)]
final_map_list[[2]] <- self_map

sgg_map <- rbindlist(final_map_list)
sgg_map <- unique(sgg_map)

cat("SGG Mapping table generated. Total unique mappings:", nrow(sgg_map), "\n")

# ----------------------------------------------------
# 매핑 함수
# ----------------------------------------------------
# ----------------------------------------------------
# 수정된 매핑 함수 (dt를 복사하여 사용)
# ----------------------------------------------------
convert_sgg <- function(dt, map_dt) {
  dt_copy <- copy(dt) # 🎯 안전한 처리를 위해 복사본 사용
    target_cols <- grep("SGG|sgg|SIGUNGU|code|CODE", names(dt_copy), value = TRUE)
  if (!("SGG_CD" %in% names(dt_copy)) && length(target_cols) > 0) {
    setnames(dt_copy, old = target_cols[1], new = "SGG_CD")
  } else if (!("SGG_CD" %in% names(dt_copy))) {
    cat("  ⚠️ Warning: SGG code column not found. Skipping SGG mapping.\n")
    return(dt) # SGG_CD가 없으면 원본 반환
  }
  
  dt_copy[, old_code := as.character(SGG_CD)]
  dt_copy[map_dt, on = .(old_code), SGG_CD := i.new_code]
  dt_copy[, old_code := NULL]
  
  dt_copy <- dt_copy[!is.na(SGG_CD)]
  return(dt_copy)
}

## =========================================
## 1. 연도별 노출 / 건강 데이터 읽기 (hourly)
## =========================================
expo_list <- list()
ha_list   <- list()

for (year_use in years) {
  cat("\n===============================\n")
  cat(" Reading year:", year_use, "\n")
  cat("===============================\n")
  
  ## -----------------------------
  ## 1-1. Exposure (total_cleaned, hourly)
  ## -----------------------------
  expo_file <- file.path(expo_dir, sprintf("%d_total_cleaned.csv", year_use))
  cat("  expo_file:", expo_file, "\n")
  if (!file.exists(expo_file)) {
    cat("  ❌ expo file not found, skip.\n")
    next
  }
  expo <- fread(expo_file)
  cat("  raw nrow(expo):", nrow(expo), "\n")
  
  # SGG 코드 매핑
  expo <- convert_sgg(expo, sgg_map)
  expo[, SGG_CD := as.character(SGG_CD)]
  
  # datetime 처리
  if (!("datetime" %in% names(expo))) {
    stop("Exposure에 datetime 컬럼이 없습니다. year=", year_use)
  }
  if (!inherits(expo$datetime, "POSIXt")) {
    expo[, datetime := ymd_hms(as.character(datetime), tz = "Asia/Seoul")]
  }
  
  expo[, date := as.Date(datetime)]
  expo[, `:=`(
    year  = year(datetime),
    month = month(datetime),
    wday  = wday(datetime),
    hour  = hour(datetime)
  )]
  expo <- expo[year == year_use]
  
  # PM 변수 선택 (PM25_wildfire_excess 우선, 없으면 PM25)
  pm_col <- intersect(c("PM25_wildfire_excess", "PM25"), names(expo))[1]
  if (is.na(pm_col)) {
    stop("PM 변수(PM25_wildfire_excess 또는 PM25)가 expo에 없습니다. year=", year_use)
  }
  cat("  ▶ Using PM column:", pm_col, "\n")
  expo[, pm := as.numeric(get(pm_col))]
  
  # temp 컬럼이 있으면 사용 (없으면 NA)
  temp_col <- intersect(c("temp", "기온(°C)", "TEMP"), names(expo))[1]
  if (!is.na(temp_col)) {
    cat("  ▶ Using TEMP column:", temp_col, "\n")
    expo[, temp := as.numeric(get(temp_col))]
  } else {
    expo[, temp := NA_real_]
  }
  
  expo_list[[as.character(year_use)]] <- expo[
    , .(SGG_CD, datetime, date, year, month, wday, hour, pm, temp)
  ]
  
  ## -----------------------------
  ## 1-2. Health (NEDIS_count_hourly)
  ## -----------------------------
  ha_file <- file.path(nedis_dir, sprintf("%d_nt_with_centroid_filled.csv", year_use))
  cat("  ha_file  :", ha_file, "\n")
  if (!file.exists(ha_file)) {
    cat("  ❌ NEDIS file not found, skip.\n")
    next
  }
  ha <- fread(ha_file)
  cat("  raw nrow(ha):", nrow(ha), "\n")
  
  ha <- convert_sgg(ha, sgg_map)
  ha[, SGG_CD := as.character(SGG_CD)]
  
  if (!("datetime" %in% names(ha))) {
    stop("NEDIS에 datetime 컬럼이 없습니다. year=", year_use)
  }
  if (!inherits(ha$datetime, "POSIXt")) {
    ha[, datetime := ymd_hms(as.character(datetime), tz = "Asia/Seoul")]
  }
  
  ha[, `:=`(
    case_datetime = datetime,
    date  = as.Date(datetime),
    year  = year(datetime),
    month = month(datetime),
    wday  = wday(datetime),
    hour  = hour(datetime)
  )]
  
  if (!("count" %in% names(ha))) {
    stop("NEDIS에 count 컬럼이 없습니다. year=", year_use)
  }
  
  # SGG×datetime별 count 합산
  ha <- ha[, .(count = sum(count, na.rm = TRUE),
               case_datetime = unique(case_datetime),
               date  = unique(date),
               year  = unique(year),
               month = unique(month),
               wday  = unique(wday),
               hour  = unique(hour)),
           by = .(SGG_CD, datetime)]
  
  ha_list[[as.character(year_use)]] <- ha[
    , .(SGG_CD, datetime, case_datetime, date, year, month, wday, hour, count)
  ]
}

## 합치기
expo_all <- rbindlist(expo_list, use.names = TRUE, fill = TRUE)
ha_all   <- rbindlist(ha_list,   use.names = TRUE, fill = TRUE)

cat("\nTotal exposure rows:", nrow(expo_all), "\n")
cat("Total health rows   :", nrow(ha_all), "\n")

## =========================================
## 2. case-crossover용 dt 만들기 (hourly, time-stratified)
## =========================================

## 사건(케이스)만: count > 0
resp <- ha_all[count > 0]
resp[, id := .I] # strata id
resp[, outcome := count]

## expo에서 필요한 컬럼만
expo_small <- expo_all[
  , .(SGG_CD, datetime, date, year, month, wday, hour, pm, temp)
]
# 🎯 [수정 1] 노출 시점 datetime의 이름을 명확하게 변경
setnames(expo_small, "datetime", "exposure_datetime")


## time-stratified 매칭:
##  - 같은 SGG_CD
##  - 같은 year, month, wday, hour (요일/월/연도/시 동일)
dt <- merge(
  resp,
  expo_small,
  by = c("SGG_CD", "year", "month", "wday", "hour"),
  all.x = TRUE,
  allow.cartesian = TRUE
)

# 🎯 [수정 2] case indicator를 새로운 exposure_datetime 기준으로 재작성
dt[, case := fifelse(case_datetime == exposure_datetime, 1L, 0L)]

setnames(dt, "count", "outcome") # 혹시 남아있을 경우 통일

## NA pm / temp 제거 (pm은 필수, temp는 있으면 같이 필터)
dt <- dt[!is.na(pm)]
# temp가 NA인 행은 그냥 써도 되지만, 완전 케이스만 쓰고 싶으면 아래 활성화
# dt <- dt[!is.na(temp)]

cat("\n nrow(dt) after filtering:", nrow(dt), "\n")
cat(" summary(pm):\n"); print(summary(dt$pm))
cat(" temp NA proportion:", mean(is.na(dt$temp)), "\n")

## =========================================
## 3. SGG별 케이스 기여도 확인 (Optional)
## =========================================
sgg_contribution <- dt[case == 1, .N, by = SGG_CD][order(-N)]
cat("\n=== 분석에 포함된 SGG 코드 기여도 (case 기준) ===\n")
cat("총 SGG_CD 수 (case > 0):", nrow(sgg_contribution), "\n")
print(head(sgg_contribution, 10))

fwrite(
  sgg_contribution,
  file.path(out_dir, "sgg_case_contribution_hourly_dlnm.csv")
)

## =========================================
## 4. DLNM cross-basis (lag 0–24시간)
## =========================================
# 시계열 순서를 맞춰주기 위해 정렬
setorder(dt, SGG_CD, case_datetime, datetime)

cb_pm <- crossbasis(
  dt$pm,
  lag   = 24,                    # 0~24시간
  argvar = list(fun = "lin"),    # PM 선형
  arglag = list(fun = "ns", df = 3)  # lag 방향 spline
)

## =========================================
## 5. clogit + DLNM (temp linear, method = "exact")
## =========================================

m_dlnm_temp <- clogit(
  case ~ cb_pm + temp + strata(id),
  data    = dt,
  weights = outcome,       # 환자 수로 가중
  method  = "exact"        # 🔥 exact likelihood
)

cat("\n=== DLNM + temp(linear), hourly, lag0–24, method = 'exact' ===\n")
print(summary(m_dlnm_temp))

## =========================================
## 6. 누적 RR (cumulative RR) 및 lag별 RR 추출
## =========================================

## (1) 누적 RR (lag 0–24, 10 μg/m³ 증가)
cp_cum <- crosspred(
  cb_pm,
  model = m_dlnm_temp,
  cen   = 0,
  at    = 10,
  cumul = TRUE
)

cRR <- cp_cum$cumfit
cat("\n=== 누적 RR (0→10 μg/m³, lag0–24h) ===\n")
print(cRR)

## (2) lag별 RR (lag 0–24, 10 μg/m³ 증가)
cp_lag <- crosspred(
  cb_pm,
  model = m_dlnm_temp,
  cen   = 0,
  at    = 10,
  cumul = FALSE
)

logRR_vec <- as.numeric(cp_lag$matfit)
lags_vec  <- 0:(length(logRR_vec) - 1)   # 0~24

lag_rr_dt <- data.table(
  lag_hour = lags_vec,
  logRR    = logRR_vec,
  RR       = exp(logRR_vec)
)

cat("\n=== lag별 RR (0–24h, at = +10 μg/m³) ===\n")
print(lag_rr_dt)

## =========================================
## 7. 결과 저장
## =========================================
out_file_lag <- file.path(out_dir, "dlnm_temp_only_hourly_lag0_24_rr.csv")
fwrite(lag_rr_dt, out_file_lag, bom = TRUE)

cat(sprintf("\n-> lag별 RR 결과가 %s 에 저장되었습니다.\n", out_file_lag))
