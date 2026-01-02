library(data.table)
library(lubridate)
library(survival)

## ==========================
## 0. 경로 설정
## ==========================
base_dir   <- "C:/Users/USER/Desktop/thesis"
expo_dir   <- file.path(base_dir, "wf_bi")
nedis_dir  <- file.path(base_dir, "NEDIS_count_daily")
out_dir    <- file.path(base_dir, "output1")
map_dir    <- file.path(base_dir, "sgg_cd") 
dir.create(out_dir, showWarnings = FALSE)

years <- 2015:2024

expo_list <- list()
ha_list   <- list()

## =========================================
## 0.5. SGG 매핑 테이블 생성 (최대 보완)
## =========================================
cat("Generating robust SGG mapping table (Target 2024 standard, all exceptions)...\n")

map_file_2022 <- file.path(map_dir, "sgg_cd_2022.csv")
map_file_2023 <- file.path(map_dir, "sgg_cd_2023.csv")
map_file_2024 <- file.path(map_dir, "sgg_cd_2024.csv")

if (!all(file.exists(map_file_2022, map_file_2023, map_file_2024))) {
  stop(paste("❌ Error: SGG mapping files not found in", map_dir))
}

# 1. 데이터 로드 및 전처리 (문자열 통일)
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

# 2. Historical (2022, 2023) 코드 통합
df_hist <- rbindlist(list(
  df_2022[, .(old_code = code, sido, sigungu)],
  df_2023[, .(old_code = code, sido, sigungu)]
), use.names = TRUE)
df_hist <- unique(df_hist)

# 3. Target List (2024년 기준)
target_list <- unique(df_2024[, .(SGG_CD = code)])

# 4. Base Mapping (sido + sigungu 이름 기준 병합)
base_map <- merge(
  df_hist,
  df_2024[, .(new_code = code, sido, sigungu)],
  by = c("sido", "sigungu"), all.x = TRUE
)

# 5. 명시적 예외 처리 (이름으로 매칭되지 않는 행정구역 변경)
# 5-1. 🎯 군위군: sido 변경 (경북 47720 -> 대구 27720)
base_map[sigungu == "군위군" & old_code == "47720" & is.na(new_code), new_code := "27720"] 

# 5-2. 🎯 인천 남구 -> 미추홀구 (이름 변경)
michuhol_code_2024 <- df_2024[sigungu %like% "미추홀", code][1]
if(!is.na(michuhol_code_2024)) {
  # 28170은 인천 남구의 옛 코드 (이름 변경 전)
  base_map[old_code == "28170" & is.na(new_code), new_code := michuhol_code_2024]
}

# 5-3. 🎯 부천시 (3구 -> 통합 -> 3구)
# 2024년 부천의 구 코드 중 하나를 대표 코드로 사용 (예: 원미구)
bucheon_target_2024 <- df_2024[sido == "경기" & sigungu %like% "원미구", code][1] 

if (!is.na(bucheon_target_2024)) {
  # [A] 2016년 이전 3개 구 코드 (41195:원미, 41197:소사, 41199:오정)
  old_3gu_codes <- c("41195", "41197", "41199")
  base_map[old_code %in% old_3gu_codes & is.na(new_code), new_code := bucheon_target_2024]
  
  # [B] 2016~2019 통합 '부천시' 코드 (df_hist에 존재할 수 있음)
  bucheon_single_codes <- df_hist[sigungu == "부천시" & sido == "경기", old_code] 
  base_map[old_code %in% bucheon_single_codes & is.na(new_code), new_code := bucheon_target_2024]
}


# 6. 최종 매핑 테이블 생성
final_map_list <- list()
final_map_list[[1]] <- base_map[!is.na(old_code) & !is.na(new_code), .(old_code, new_code)]
existing_olds <- final_map_list[[1]]$old_code
self_map <- target_list[!SGG_CD %in% existing_olds, .(old_code = SGG_CD, new_code = SGG_CD)]
final_map_list[[2]] <- self_map

sgg_map <- rbindlist(final_map_list)
sgg_map <- unique(sgg_map)

cat("SGG Mapping table generated. Total unique mappings:", nrow(sgg_map), "\n")


# ----------------------------------------------------
# 매핑 함수 (SGG_CD를 명시적으로 문자열로 사용)
# ----------------------------------------------------
convert_sgg <- function(dt, map_dt) {
  target_cols <- grep("SGG|sgg|SIGUNGU|code|CODE", names(dt), value = TRUE)
  if (!("SGG_CD" %in% names(dt)) && length(target_cols) > 0) {
    setnames(dt, old = target_cols[1], new = "SGG_CD")
  } else if (!("SGG_CD" %in% names(dt))) {
    cat("  ⚠️ Warning: SGG code column not found. Skipping SGG mapping.\n")
    return(dt)
  }
  
  dt[, old_code := as.character(SGG_CD)]
  dt[map_dt, on = .(old_code), SGG_CD := i.new_code]
  dt[, old_code := NULL] 
  
  dt <- dt[!is.na(SGG_CD)]
  
  return(dt)
}

## ==================================================
## 1. 데이터 읽기 (SGG 매핑 적용)
## ==================================================

## ----------------------------------------
## 1-1. Exposure (wf_bi) 데이터 읽기
## ----------------------------------------
for (year_use in years) {
  cat("\n===============================\n")
  cat(" Reading Exposure year:", year_use, "\n")
  cat("===============================\n")
  
  # 🎯 FIX: Exposure 파일명 수정
  expo_file <- file.path(expo_dir, paste0(year_use, "_wfbi_with_sgg.csv"))
  if (!file.exists(expo_file)) { cat("  ❌ Error: Exposure file not found for year", year_use, "\n"); next }
  
  expo <- fread(expo_file)
  if (nrow(expo) == 0) { cat("  ❌ Error: Exposure data is empty for year", year_use, "\n"); next }
  
  if (!is.character(expo$SGG_CD)) { expo[, SGG_CD := as.character(SGG_CD)] }
  expo <- convert_sgg(expo, sgg_map)
  
  ## ---- 날짜 및 변수 설정 ----
  expo[, date := as.Date(date, format = "%Y-%m-%d")] 
  expo <- expo[year(date) == year_use]
  expo[, wf := smoke_day]
  expo[, year := year_use]
  
  expo_list[[as.character(year_use)]] <- expo[, .(date, SGG_CD, wf, year)]
}

## ----------------------------------------
## 1-2. Health (NEDIS_count_daily) 데이터 읽기
## ----------------------------------------
for (year_use in years) {
  cat("\n===============================\n")
  cat(" Reading Health year:", year_use, "\n")
  cat("===============================\n")
  
  ha_file <- file.path(nedis_dir, paste0(year_use, "_daily_with_centroid.csv"))
  if (!file.exists(ha_file)) { cat("  ❌ Error: Health file not found for year", year_use, "\n"); next }
  
  ha <- NULL 
  tryCatch({
    ha <- fread(ha_file, colClasses = c(SGG_CD="character"))
  }, error = function(e) {
    cat("  ❌ Error reading Health file for year", year_use, ":", conditionMessage(e), "\n")
  })
  
  if (is.null(ha) || nrow(ha) == 0) { cat("  ❌ Error: Health data (ha) is empty or failed to load for year", year_use, "\n"); next }
  
  if ("date" %in% names(ha) && "daily_count" %in% names(ha)) {
    setnames(ha, old = "date", new = "case_date")
    setnames(ha, old = "daily_count", new = "respiratory")
  } else {
    cat("  ❌ Error: Required columns (date or daily_count) not found in Health data for year", year_use, "\n"); next
  }
  
  ha <- convert_sgg(ha, sgg_map)
  if (!is.character(ha$SGG_CD)) { ha[, SGG_CD := as.character(SGG_CD)] }
  
  ha[, case_date := as.Date(case_date, format = "%Y-%m-%d")]
  
  ha <- ha[!is.na(case_date) & year(case_date) == year_use]
  ha[, year := year_use]
  
  ha[, `:=`( wday  = wday(case_date), month = month(case_date) )]
  
  ha_list[[as.character(year_use)]] <- ha[, .(SGG_CD, respiratory, case_date, year, wday, month)]
}

## ==========================
## 2. 전체 연도 합치기
## ==========================
expo_all <- rbindlist(expo_list, use.names = TRUE, fill = TRUE)
ha_all   <- rbindlist(ha_list,   use.names = TRUE, fill = TRUE)

cat("\nTotal exposure rows:", nrow(expo_all), "\n")
cat("Total health rows   :", nrow(ha_all), "\n")

if (nrow(expo_all) == 0 || nrow(ha_all) == 0) {
  stop("Error: Total exposure or health data is empty after loading.")
}

## ==========================
## 3. Pooled case-crossover (RR, PAN, PAF 계산 - Lag 7까지)
## ==========================
lags <- c("same_day", paste0("lag", 1:7)) 
pooled_results <- list() 

for (lg in lags) {
  cat("\n-----------------------------\n")
  cat("Running pooled lag:", lg, "\n")
  cat("-----------------------------\n")
  
  ## 3-1. exposure + lag
  new <- copy(expo_all)[order(SGG_CD, date)]
  if (lg != "same_day") {
    nn <- as.numeric(gsub("lag","",lg)) 
    new[, wf := shift(wf, n = nn, type = "lag"), by = SGG_CD]
  }
  
  new[, `:=`(
    wday  = wday(date),
    month = month(date),
    year  = year(date)
  )]
  
  ## 3-2. resp_all & id 
  resp_all <- ha_all[respiratory > 0]
  resp_all[, id := .I]
  
  ## 3-3. merge → dt_all (Matching)
  dt_all <- merge(
    resp_all,
    new,
    by = c("SGG_CD","wday","month","year"),
    all.x = TRUE,
    allow.cartesian = TRUE
  )
  
  dt_all[, case := fifelse(case_date == date, 1L, 0L)]
  dt_all <- dt_all[!is.na(wf)] 
  setnames(dt_all, "respiratory", "outcome") 
  
  cat("  nrow(dt_all) =", nrow(dt_all),
      ", sum(case) =", sum(dt_all$case), "\n")
  
  if (nrow(dt_all) == 0 || sum(dt_all$case) == 0) {
    cat("  ❌ Skipping lag", lg, ": No valid case-control pairs found.\n")
    next 
  }
  
  ## 3-4. clogit (RR, PAN, PAF 계산)
  fit <- clogit(
    case ~ wf + strata(id), 
    data    = dt_all,
    weights = outcome,
    method  = "approximate"
  )
  
  summ <- summary(fit)
  ci   <- summ$conf.int[1, c(1,3,4)]
  RR    <- ci[1]
  RR_ll <- ci[2]
  RR_ul <- ci[3]
  
  exposed_cases <- dt_all[case == 1 & wf == 1, sum(outcome, na.rm = TRUE)]
  
  paf    <- 1 - (1/RR)
  paf_ll <- 1 - (1/RR_ll)
  paf_ul <- 1 - (1/RR_ul)
  
  pan    <- exposed_cases * paf
  pan_ll <- exposed_cases * paf_ll
  pan_ul <- exposed_cases * paf_ul
  
  pooled_results[[lg]] <- data.table(
    lag     = lg,
    RR      = RR,
    RR_ll   = RR_ll,
    RR_ul   = RR_ul,
    exposed_n = exposed_cases,
    PAF     = paf,
    PAN     = pan,
    PAN_ll  = pan_ll,
    PAN_ul  = pan_ul
  )
}

## =========================================
## 4. 결과 저장
## =========================================
final_pooled <- rbindlist(pooled_results)

print(final_pooled)

fwrite(final_pooled, 
       file.path(out_dir, "Pooled_Statewide_RR_PAN_2015_2024.csv"), 
       bom = TRUE)

cat("\nSaved to:", file.path(out_dir, "Pooled_Statewide_RR_PAN_2015_2024.csv"), "\n")