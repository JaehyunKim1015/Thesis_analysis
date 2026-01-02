library(data.table)
library(lubridate)
library(survival)

## ==========================================
## 0. 경로 및 연도 설정
## ==========================================
expo_dir  <- "C:/Users/USER/Desktop/thesis/wildfire/total_clean"
nedis_dir <- "C:/Users/USER/Desktop/thesis/NEDIS_count_hourly"
out_dir   <- "C:/Users/USER/Desktop/thesis/output3"
map_dir   <- "C:/Users/USER/Desktop/thesis/sgg_cd"

dir.create(out_dir, showWarnings = FALSE)

years <- 2015:2024   # 필요시 조정

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

# 3. Target List (2024년 기준, 코드만)
target_list <- unique(df_2024[, .(SGG_CD = code)])

# 4. Base Mapping (sido + sigungu 이름 기준 병합)
base_map <- merge(
  df_hist,
  df_2024[, .(new_code = code, sido, sigungu)],
  by = c("sido", "sigungu"),
  all.x = TRUE
)

# 5. 명시적 예외 처리
# 5-1. 군위군: 경북(47720) -> 대구(27720)
base_map[sigungu == "군위군" & old_code == "47720" & is.na(new_code), new_code := "27720"] 

# 5-2. 인천 남구 -> 미추홀구
michuhol_code_2024 <- df_2024[sigungu %like% "미추홀", code][1]
if (!is.na(michuhol_code_2024)) {
  base_map[old_code == "28170" & is.na(new_code), new_code := michuhol_code_2024]
}

# 5-3. 부천시 (3구 → 통합 → 3구) 대표코드(예: 원미구)
bucheon_target_2024 <- df_2024[sido == "경기" & sigungu %like% "원미구", code][1]

if (!is.na(bucheon_target_2024)) {
  old_3gu_codes <- c("41195", "41197", "41199")
  base_map[old_code %in% old_3gu_codes & is.na(new_code), new_code := bucheon_target_2024]
  
  bucheon_single_codes <- df_hist[sigungu == "부천시" & sido == "경기", old_code]
  base_map[old_code %in% bucheon_single_codes & is.na(new_code), new_code := bucheon_target_2024]
}

# 6. 최종 매핑 테이블 생성
final_map_list <- list()
final_map_list[[1]] <- base_map[!is.na(old_code) & !is.na(new_code),
                                .(old_code, new_code)]
existing_olds <- final_map_list[[1]]$old_code
self_map <- target_list[!SGG_CD %in% existing_olds,
                        .(old_code = SGG_CD, new_code = SGG_CD)]
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

## ==========================================
## 1. 한 연도씩 case-crossover용 데이터 생성 (hourly)
##    - wildfire_affected 기준
## ==========================================
make_cc_data_one_year <- function(YR) {
  cat("\n=====================================\n")
  cat("      ", YR, "년 데이터 준비 시작\n")
  cat("=====================================\n")
  
  ## ---------- 1) 노출 데이터 ----------
  exp_path <- file.path(expo_dir, sprintf("%d_total_cleaned.csv", YR))
  if (!file.exists(exp_path)) {
    cat("  ❌ Exposure 파일 없음:", exp_path, "\n")
    return(NULL)
  }
  
  exp_dt <- fread(exp_path)
  cat("  ▶ EXP 컬럼:", paste(names(exp_dt), collapse = ", "), "\n")
  
  # ✅ SGG 매핑
  exp_dt <- convert_sgg(exp_dt, sgg_map)
  
  # 형식 정리
  exp_dt[, SGG_CD := as.character(SGG_CD)]
  
  if (!("datetime" %in% names(exp_dt))) {
    stop("Exposure에 datetime 컬럼이 없습니다. (year=", YR, ")")
  }
  
  # datetime을 POSIXct로 (이미 POSIX면 그대로)
  if (!inherits(exp_dt$datetime, "POSIXt")) {
    exp_dt[, datetime := ymd_hms(as.character(datetime), tz = "Asia/Seoul")]
  }
  
  cat("  ▶ EXP datetime NA 개수:", sum(is.na(exp_dt$datetime)), "\n")
  exp_dt <- exp_dt[!is.na(datetime)]
  
  # wildfire_affected 0/1 보장
  if (!("wildfire_affected" %in% names(exp_dt))) {
    stop("Exposure에 wildfire_affected 컬럼이 없습니다. (year=", YR, ")")
  }
  exp_dt[, wildfire_affected := as.integer(wildfire_affected > 0)]
  
  # 시간 파생변수
  exp_dt[, `:=`(
    year  = year(datetime),
    month = month(datetime),
    dow   = wday(datetime),
    hour  = hour(datetime)
  )]
  
  exp_dt <- exp_dt[year == YR]
  cat("  ▶ EXP 행수 (연도 필터 후):", nrow(exp_dt), "\n")
  
  ## ---------- 2) NEDIS 데이터 ----------
  ed_path <- file.path(nedis_dir, sprintf("%d_nt_with_centroid_filled.csv", YR))
  if (!file.exists(ed_path)) {
    cat("  ❌ NEDIS 파일 없음:", ed_path, "\n")
    return(NULL)
  }
  
  ed_dt <- fread(ed_path)
  cat("  ▶ NEDIS 컬럼:", paste(names(ed_dt), collapse = ", "), "\n")
  
  # ✅ SGG 매핑
  ed_dt <- convert_sgg(ed_dt, sgg_map)
  
  ed_dt[, SGG_CD := as.character(SGG_CD)]
  
  if (!("datetime" %in% names(ed_dt))) {
    stop("NEDIS에 datetime 컬럼이 없습니다. (year=", YR, ")")
  }
  
  if (!inherits(ed_dt$datetime, "POSIXt")) {
    ed_dt[, datetime := ymd_hms(as.character(datetime), tz = "Asia/Seoul")]
  }
  
  cat("  ▶ NEDIS datetime NA 개수:", sum(is.na(ed_dt$datetime)), "\n")
  ed_dt <- ed_dt[!is.na(datetime)]
  
  # 해당 연도만
  ed_dt[, year := year(datetime)]
  ed_dt <- ed_dt[year == YR]
  
  if (!("count" %in% names(ed_dt))) {
    stop("NEDIS에 count 컬럼이 없습니다. (year=", YR, ")")
  }
  
  # SGG_CD + datetime별 count 합산
  ed_dt <- ed_dt[, .(count = sum(count, na.rm = TRUE)), by = .(SGG_CD, datetime)]
  cat("  ▶ NEDIS 정리 후 행수:", nrow(ed_dt), "\n")
  
  ## ---------- 3) 노출 + NEDIS 병합 ----------
  merged <- merge(
    exp_dt[, .(SGG_CD, datetime, wildfire_affected, year, month, dow, hour)],
    ed_dt,
    by = c("SGG_CD", "datetime"),
    all.x = TRUE
  )
  
  merged[is.na(count), count := 0L]
  
  cat("  ▶ 병합 후 행수:", nrow(merged), "\n")
  cat("  ▶ count > 0 (케이스) 행수:", sum(merged$count > 0), "\n")
  
  ## ---------- 4) case-crossover 변수 ----------
  merged[, `:=`(
    case    = as.integer(count > 0),
    weight  = ifelse(count > 0, count, 1L),
    stratum = interaction(SGG_CD, year, month, dow, hour, drop = TRUE)
  )]
  
  cc_dt <- merged[!is.na(wildfire_affected)]
  
  cat("  ▶ 최종 분석용 행수:", nrow(cc_dt), "\n")
  cat("  ▶ 케이스 수:", sum(cc_dt$case == 1), "\n")
  
  return(cc_dt)
}

## ==========================================
## 2. 연도별 데이터 생성 & 합치기
## ==========================================
cc_list <- lapply(years, make_cc_data_one_year)
cc_list <- cc_list[!sapply(cc_list, is.null)]

if (length(cc_list) == 0) {
  stop("❌ 유효한 연도 데이터가 없습니다.")
}

cc_all <- rbindlist(cc_list, use.names = TRUE, fill = TRUE)

cat("\n===== 전체 기간 데이터 행수:", nrow(cc_all), "=====\n")
cat("===== 전체 케이스 수:", sum(cc_all$case == 1), "=====\n")

cc_all[, wildfire_affected := as.integer(wildfire_affected > 0)]

## =========================================
## 3. 분석 모델링 (Lag 0, 1, 12 / SGG별, Time-stratified C-C)
##    - 노출: wildfire_affected
##    - 방법: exact clogit
## =========================================
lags <- c(0, 1, 12) # 분석할 시간 lag (hours)
res_lag_list <- list()

# cc_all 데이터테이블 복사 및 정렬 (Lag 계산을 위해)
cc_ordered <- copy(cc_all)[order(SGG_CD, datetime)]

# 실제 데이터에 존재하는 SGG 코드만 사용
all_sgg_codes <- unique(cc_ordered$SGG_CD)
cat("총 SGG 개수:", length(all_sgg_codes), "\n")
cat("총 케이스 수:", sum(cc_ordered$case), "\n")

for (lg in lags) {
  cat("\n============================\n")
  cat("  ▶ Running Lag (hours):", lg, "\n")
  cat("============================\n")
  
  # 3-1. Lag 변수 생성 (wildfire_affected 기준)
  if (lg == 0) {
    cc_ordered[, wf_lag := wildfire_affected]
  } else {
    cc_ordered[, wf_lag := shift(wildfire_affected, n = lg, type = "lag"), by = SGG_CD]
  }
  
  # Lag 변수 결측치 제거
  dt_lag <- cc_ordered[!is.na(wf_lag)]
  
  # 3-2. SGG별 clogit (Time-stratified Case-Crossover)
  res_sgg <- vector("list", length(all_sgg_codes))
  
  for (i in seq_along(all_sgg_codes)) {
    sgg <- all_sgg_codes[i]
    
    # 해당 시군구에서 case가 있었던 stratum만 추출
    foo <- dt_lag[SGG_CD == sgg & case == 1]
    if (nrow(foo) == 0) {
      res_sgg[[i]] <- data.table(
        SGG_CD   = sgg,
        lag_hour = lg,
        RR       = NA_real_,
        RR_ll    = NA_real_,
        RR_ul    = NA_real_,
        n_case   = 0L,
        note     = "No Case"
      )
      next
    }
    
    unique_strata <- unique(foo$stratum)
    
    # 해당 SGG + strata에 포함되는 모든 시간 (case + control)
    dt_sgg_cc <- dt_lag[SGG_CD == sgg & stratum %in% unique_strata]
    
    n_case_val <- sum(dt_sgg_cc$case)
    
    # clogit 실행 조건:
    # - 케이스 존재
    # - 노출변수(wf_lag)가 0/1 양쪽 다 존재
    # - strata 최소 1개 이상
    if (n_case_val > 0 &&
        length(unique(dt_sgg_cc$wf_lag)) >= 2 &&
        length(unique(dt_sgg_cc$stratum)) > 0) {
      
      fit <- try(
        clogit(
          case ~ wf_lag + strata(stratum),
          data    = dt_sgg_cc,
          weights = weight,      # 환자 수(응급실 방문 수)로 가중
          method  = "exact"      # exact 방법
        ),
        silent = TRUE
      )
      
      if (!inherits(fit, "try-error")) {
        ci <- summary(fit)$conf.int[1, c("exp(coef)", "lower .95", "upper .95")]
        res_sgg[[i]] <- data.table(
          SGG_CD   = sgg,
          lag_hour = lg,
          RR       = ci[1],
          RR_ll    = ci[2],
          RR_ul    = ci[3],
          n_case   = n_case_val,
          note     = "Success"
        )
      } else {
        res_sgg[[i]] <- data.table(
          SGG_CD   = sgg,
          lag_hour = lg,
          RR       = NA_real_,
          RR_ll    = NA_real_,
          RR_ul    = NA_real_,
          n_case   = n_case_val,
          note     = "Model Error"
        )
      }
      
    } else {
      res_sgg[[i]] <- data.table(
        SGG_CD   = sgg,
        lag_hour = lg,
        RR       = NA_real_,
        RR_ll    = NA_real_,
        RR_ul    = NA_real_,
        n_case   = n_case_val,
        note     = "Insufficient Data"
      )
    }
    
    if (i %% 50 == 0) cat(".")
  }
  
  res_lag_list[[as.character(lg)]] <- rbindlist(res_sgg)
}

# 최종 결과 합치기
res_all <- rbindlist(res_lag_list)

## =========================================
## 4. SGG 이름 붙여서 저장
## =========================================
# 2024 SGG 이름 정보 다시 로드
map_file_2024 <- file.path(map_dir, "sgg_cd_2024.csv")
map_df_2024 <- fread(map_file_2024)
setnames(map_df_2024, c("code", "sido", "sigungu"))
map_df_2024[, code := as.character(code)]
map_df_2024[, `:=`(sido = remove_space(sido), sigungu = remove_space(sigungu))]
target_sgg_names <- unique(map_df_2024[, .(
  SGG_CD  = code,
  SIDO_NM = sido,
  SGG_NM  = sigungu
)])

final_res <- merge(
  target_sgg_names,  # SGG_CD, SIDO_NM, SGG_NM
  res_all,
  by = "SGG_CD",
  all.x = TRUE
)

final_res <- final_res[order(SIDO_NM, SGG_NM, lag_hour)]

cat("\nTotal Rows:", nrow(final_res), "\n")
cat("성공적으로 모델링된 SGG 개수:",
    sum(final_res$note == "Success", na.rm = TRUE), "\n")
print(final_res[note == "Success" & lag_hour == 0][1:5])

out_file <- file.path(out_dir, "SGG_level_RR_2015_2024_wildfire_affected_cc_hourly.csv")
fwrite(final_res, out_file, bom = TRUE)
cat("\nSaved to:", out_file, "\n")
