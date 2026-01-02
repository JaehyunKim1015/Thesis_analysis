library(data.table)
library(lubridate)
library(survival)

## ==========================================
## 0. 경로 및 연도 설정
## ==========================================
base_dir  <- "C:/Users/USER/Desktop/thesis"

expo_dir  <- file.path(base_dir, "wildfire", "total_clean")      # ✅ total_clean 사용
nedis_dir <- file.path(base_dir, "NEDIS_count_hourly")
out_dir   <- file.path(base_dir, "output3")
map_dir   <- file.path(base_dir, "sgg_cd")

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

# 3. Target List (2024년 기준)
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

## ==========================================
## 1. 한 연도씩 case-crossover용 데이터 생성 (시간단위)
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
  
  # ✅ SGG 매핑 먼저 적용
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
  
  # ✅ SGG 매핑 적용
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

## ==========================================
## 3. lag 0 (동시간) time-stratified case-crossover + PAN
## ==========================================
model_wf <- clogit(
  case ~ wildfire_affected + strata(stratum),
  weights = weight,
  data    = cc_all
  )


print(summary(model_wf))

beta <- coef(model_wf)["wildfire_affected"]
ci   <- confint(model_wf)["wildfire_affected", ]

OR  <- exp(beta)
LCL <- exp(ci[1])
UCL <- exp(ci[2])

cat("\n=== Wildfire_affected 효과 (lag 0, 단변량) ===\n")
cat("OR :", OR, "\n")
cat("95% CI:", LCL, "~", UCL, "\n")

## ---------- PAN 계산 ----------
total_cases   <- sum(cc_all$case == 1, na.rm = TRUE)
exposed_cases <- sum(cc_all$case == 1 & cc_all$wildfire_affected == 1, na.rm = TRUE)
Pe            <- ifelse(total_cases > 0, exposed_cases / total_cases, NA_real_)

PAF   <- Pe * (OR  - 1) / OR
PAF_L <- Pe * (LCL - 1) / LCL
PAF_U <- Pe * (UCL - 1) / UCL

PAN   <- PAF   * total_cases
PAN_L <- PAF_L * total_cases
PAN_U <- PAF_U * total_cases

cat("\n=== PAN (wildfire_affected, lag 0) ===\n")
cat("총 케이스 수:", total_cases, "\n")
cat("노출 케이스 수:", exposed_cases, "\n")
cat("Pe (케이스 중 노출 비율):", Pe, "\n")
cat("PAF:", PAF, " (95% CI:", PAF_L, "~", PAF_U, ")\n")
cat("PAN (기여사례수):", PAN, " (95% CI:", PAN_L, "~", PAN_U, ")\n")

res_main <- data.table(
  exposure      = "wildfire_affected",
  lag_hours     = 0,
  OR            = OR,
  LCL           = LCL,
  UCL           = UCL,
  total_cases   = total_cases,
  exposed_cases = exposed_cases,
  Pe            = Pe,
  PAF           = PAF,
  PAF_L         = PAF_L,
  PAF_U         = PAF_U,
  PAN           = PAN,
  PAN_L         = PAN_L,
  PAN_U         = PAN_U
)

fwrite(res_main,
       file.path(out_dir, "cc_wildfire_affected_lag0_result_with_PAN.csv"),
       bom = TRUE)

## ==========================================
## 4. lag 0~12시간 OR + PAN 계산
## ==========================================
max_lag <- 12
setorder(cc_all, SGG_CD, datetime)

# lag 변수 생성 (wildfire_affected 기준)
for (k in 0:max_lag) {
  vname <- paste0("wf_lag", k)
  if (k == 0) {
    cc_all[, (vname) := wildfire_affected]
  } else {
    cc_all[, (vname) := shift(wildfire_affected, n = k, type = "lag"), by = SGG_CD]
  }
}

lag_results <- list()

for (k in 0:max_lag) {
  cat("\n-----------------------------\n")
  cat("  ▶ lag (hours) =", k, "\n")
  cat("-----------------------------\n")
  
  vname <- paste0("wf_lag", k)
  dt <- cc_all[!is.na(get(vname))]
  if (nrow(dt) == 0) {
    cat("  ❌ 데이터 없음, 건너뜀\n")
    next
  }
  
  dt[, wf := get(vname)]
  
  tab_global <- table(dt$wf, dt$case)
  cat("  ▶ wf x case 분포:\n")
  print(tab_global)
  
  total_cases_k   <- sum(dt$case == 1, na.rm = TRUE)
  exposed_cases_k <- sum(dt$case == 1 & dt$wf == 1, na.rm = TRUE)
  
  # 변이 부족한 경우
  if (nrow(tab_global) < 2 || ncol(tab_global) < 2) {
    cat("  ❌ 변이 부족 → OR/PAN 추정 불가, NA로 기록\n")
    lag_results[[as.character(k)]] <- data.table(
      lag_hours     = k,
      OR            = NA_real_,
      LCL           = NA_real_,
      UCL           = NA_real_,
      exposed_n     = exposed_cases_k,
      total_cases   = total_cases_k,
      exposed_cases = exposed_cases_k,
      Pe            = NA_real_,
      PAF           = NA_real_,
      PAF_L         = NA_real_,
      PAF_U         = NA_real_,
      PAN           = NA_real_,
      PAN_L         = NA_real_,
      PAN_U         = NA_real_
    )
    next
  }
  
  fit <- try(
    clogit(
      case ~ wf + strata(stratum),
      weights = weight,
      data    = dt
    ),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) {
    cat("  ❌ clogit 실패:", fit, "\n")
    lag_results[[as.character(k)]] <- data.table(
      lag_hours     = k,
      OR            = NA_real_,
      LCL           = NA_real_,
      UCL           = NA_real_,
      exposed_n     = exposed_cases_k,
      total_cases   = total_cases_k,
      exposed_cases = exposed_cases_k,
      Pe            = NA_real_,
      PAF           = NA_real_,
      PAF_L         = NA_real_,
      PAF_U         = NA_real_,
      PAN           = NA_real_,
      PAN_L         = NA_real_,
      PAN_U         = NA_real_
    )
    next
  }
  
  summ   <- summary(fit)
  ci_row <- summ$conf.int[1, ]
  
  OR_k  <- ci_row["exp(coef)"]
  LCL_k <- ci_row["lower .95"]
  UCL_k <- ci_row["upper .95"]
  
  Pe_k <- ifelse(total_cases_k > 0, exposed_cases_k / total_cases_k, NA_real_)
  
  PAF_k   <- Pe_k * (OR_k  - 1) / OR_k
  PAF_L_k <- Pe_k * (LCL_k - 1) / LCL_k
  PAF_U_k <- Pe_k * (UCL_k - 1) / UCL_k
  
  PAN_k   <- PAF_k   * total_cases_k
  PAN_L_k <- PAF_L_k * total_cases_k
  PAN_U_k <- PAF_U_k * total_cases_k
  
  lag_results[[as.character(k)]] <- data.table(
    lag_hours     = k,
    OR            = OR_k,
    LCL           = LCL_k,
    UCL           = UCL_k,
    exposed_n     = exposed_cases_k,
    total_cases   = total_cases_k,
    exposed_cases = exposed_cases_k,
    Pe            = Pe_k,
    PAF           = PAF_k,
    PAF_L         = PAF_L_k,
    PAF_U         = PAF_U_k,
    PAN           = PAN_k,
    PAN_L         = PAN_L_k,
    PAN_U         = PAN_U_k
  )
}

lag_results_dt <- rbindlist(lag_results, use.names = TRUE, fill = TRUE)[order(lag_hours)]
print(lag_results_dt)

fwrite(lag_results_dt,
       file.path(out_dir, "cc_wildfire_affected_lag0_12_results_with_PAN.csv"),
       bom = TRUE)

cat("\n✅ lag 0~12시간 wildfire_affected OR + PAN 결과 저장 완료\n")
