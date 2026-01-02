library(data.table)
library(lubridate)
library(survival)

## =========================================
## 0. 경로 설정
## =========================================
base_dir   <- "C:/Users/USER/Desktop/thesis"
# SGG 매핑 파일 경로 (템플릿 유지)
map_dir    <- file.path(base_dir, "sgg_cd")
expo_dir   <- file.path(base_dir, "wf_bi")
# Health 데이터 경로 수정 (daily_count 파일명에 맞춤)
nedis_dir  <- file.path(base_dir, "NEDIS_count_daily")
out_dir    <- file.path(base_dir, "output1")

if(!dir.exists(out_dir)) dir.create(out_dir)
# 🎯 분석 기간 2024년까지 확장
years <- 2015:2024

expo_list <- list()
ha_list   <- list()

## =========================================
## 1. 강력한 매핑 테이블 생성 (2024년 표준)
## =========================================
cat("Generating robust SGG mapping table (Target 2024 standard)...\n")

map_file_2022 <- file.path(map_dir, "sgg_cd_2022.csv")
map_file_2023 <- file.path(map_dir, "sgg_cd_2023.csv")
map_file_2024 <- file.path(map_dir, "sgg_cd_2024.csv")

if (!all(file.exists(map_file_2022, map_file_2023, map_file_2024))) {
  stop(paste("❌ Error: SGG mapping files not found in", map_dir, ". Check the map_dir variable and file names."))
}

df_2022 <- fread(map_file_2022)
df_2023 <- fread(map_file_2023)
df_2024 <- fread(map_file_2024)

setnames(df_2022, c("code", "sido", "sigungu"))
setnames(df_2023, c("code", "sido", "sigungu"))
setnames(df_2024, c("code", "sido", "sigungu"))

# 🎯 FIX: 모든 코드를 문자열로 통일하여 타입 불일치 문제 해결
df_2022[, code := as.character(code)]
df_2023[, code := as.character(code)]
df_2024[, code := as.character(code)]

# [핵심 1] 공백 제거 함수 (띄어쓰기 오류 방지)
remove_space <- function(x) gsub("\\s+", "", x)
df_2022[, `:=`(sido = remove_space(sido), sigungu = remove_space(sigungu))]
df_2023[, `:=`(sido = remove_space(sido), sigungu = remove_space(sigungu))]
df_2024[, `:=`(sido = remove_space(sido), sigungu = remove_space(sigungu))]

# Historical Data 통합 (2022 & 2023)
df_hist <- rbindlist(list(
  df_2022[, .(old_code = code, sido, sigungu)],
  df_2023[, .(old_code = code, sido, sigungu)]
), use.names = TRUE)
df_hist <- unique(df_hist)

# 1-1. Target List (2024년 기준) 확정
target_list <- unique(df_2024[, .(SGG_CD = code, SIDO_NM = sido, SGG_NM = sigungu)])
target_list[, SGG_CD := as.character(SGG_CD)]

# 1-2. 기본 매핑 (이름 기준: Historical -> 2024)
base_map <- merge(
  df_hist, 
  df_2024[, .(new_code = code, sido, sigungu)], 
  by = c("sido", "sigungu"), all.x = TRUE
)

# 1-3. 필수 예외 처리 (군위, 인천, 부천)
base_map[old_code == "47720" & is.na(new_code), new_code := "27720"] # 군위군

# 인천 남구 -> 미추홀구
michuhol_code_2024 <- df_2024[sigungu %like% "미추홀", code][1]
if(!is.na(michuhol_code_2024)) base_map[old_code == "28170" & is.na(new_code), new_code := michuhol_code_2024]

# 부천시 통합/분할
bucheon_representative_2024 <- df_2024[sido == "경기" & sigungu %like% "부천시", code][1]
bucheon_old_codes <- unique(c(
  df_hist[sigungu == "부천시", old_code], 
  df_2023[sigungu %in% c("원미구", "소사구", "오정구"), code]
))
bucheon_old_codes <- bucheon_old_codes[!bucheon_old_codes %in% df_2024[sigungu %like% "부천", code]] 

if(!is.na(bucheon_representative_2024)) {
  base_map[old_code %in% bucheon_old_codes & is.na(new_code), new_code := bucheon_representative_2024]
}

# 1-4. [핵심 2] 변경되지 않은 코드들도 누락되지 않도록 '자기 자신' 매핑 추가
final_map_list <- list()
final_map_list[[1]] <- base_map[!is.na(old_code) & !is.na(new_code), .(old_code, new_code)]

existing_mapped_olds <- final_map_list[[1]]$old_code
self_map <- target_list[!SGG_CD %in% existing_mapped_olds, .(old_code = SGG_CD, new_code = SGG_CD)]
final_map_list[[2]] <- self_map

sgg_map <- rbindlist(final_map_list)
sgg_map <- unique(sgg_map)

cat("SGG Mapping table generated. Total unique mappings:", nrow(sgg_map), "\n")


# 매핑 적용 함수
convert_sgg <- function(dt, map_dt) {
  if (!is.character(dt$SGG_CD)) { dt[, SGG_CD := as.character(SGG_CD)] }
  
  dt[, old_code := SGG_CD]
  dt[map_dt, on = .(old_code), SGG_CD := i.new_code]
  dt[, old_code := NULL] 
  
  dt <- dt[!is.na(SGG_CD)]
  return(dt)
}

# 안전한 컬럼명 변경 함수 (템플릿의 함수 유지)
safe_rename_sgg <- function(dt) {
  target_cols <- grep("SGG|sgg|SIGUNGU|code|CODE", names(dt), value = TRUE)
  if ("SGG_CD" %in% names(dt)) {
    return(dt)
  } else if (length(target_cols) > 0) {
    setnames(dt, old = target_cols[1], new = "SGG_CD")
  } else {
    stop("파일에 SGG(지역코드) 관련 컬럼이 없습니다. 확인해주세요.")
  }
  return(dt)
}


## =========================================
## 2. 데이터 로드 및 집계
## =========================================
for (year_use in years) {
  cat("Processing year:", year_use, "...")
  
  ## -----------------------------
  ## 2-1. Exposure 데이터 처리
  ## -----------------------------
  expo_file <- file.path(expo_dir, paste0(year_use, "_wfbi_with_sgg.csv"))
  if (!file.exists(expo_file)) { cat(" ❌ Expo file not found.\n"); next }
  
  expo <- fread(expo_file)
  
  expo <- safe_rename_sgg(expo)
  
  wf_col <- grep("smoke_day|wf", names(expo), ignore.case = TRUE, value = TRUE)[1]
  if(!is.na(wf_col)) setnames(expo, old = wf_col, new = "wf")
  
  expo[, date := as.Date(date, format = "%Y-%m-%d")]
  expo <- expo[year(date) == year_use]
  
  expo <- convert_sgg(expo, sgg_map) 
  
  expo_agg <- expo[, .(PM25 = mean(PM25, na.rm=TRUE), wf = max(wf, na.rm=TRUE)), by=.(SGG_CD, date)]
  expo_agg[, year := year_use]
  expo_list[[as.character(year_use)]] <- expo_agg
  
  ## -----------------------------
  ## 2-2. Health 데이터 처리
  ## -----------------------------
  # 🎯 FIX: Health 파일명/경로를 _daily_with_centroid.csv 파일로 수정
  ha_file <- file.path(nedis_dir, paste0(year_use, "_daily_with_centroid.csv"))
  if (!file.exists(ha_file)) { cat(" ❌ Health file not found.\n"); next }
  
  # SGG_CD가 문자열로 로드되도록 설정
  ha <- fread(ha_file, colClasses = c(SGG_CD="character"))
  
  # 컬럼 이름 변경 (daily_count -> respiratory, date -> case_date)
  if ("daily_count" %in% names(ha) && "date" %in% names(ha)) {
    setnames(ha, old = "daily_count", new = "respiratory")
    setnames(ha, old = "date", new = "case_date")
  } else {
    cat(" ❌ Health file columns (daily_count/date) missing.\n"); next
  }
  
  # 날짜 포맷 변환 (2024_daily_with_centroid.csv 파일 형식 가정)
  ha[, case_date := as.Date(case_date, format = "%Y-%m-%d")]
  
  ha <- ha[year(case_date) == year_use]
  
  ha <- convert_sgg(ha, sgg_map)
  
  ha_agg <- ha[, .(respiratory = sum(respiratory, na.rm=TRUE)), by=.(SGG_CD, case_date)]
  ha_agg[, `:=`(year=year_use, wday=wday(case_date), month=month(case_date))]
  
  ha_list[[as.character(year_use)]] <- ha_agg
  cat(" Done.\n")
}

cat("All data loaded successfully.\n")
expo_all <- rbindlist(expo_list, fill=TRUE)
ha_all   <- rbindlist(ha_list, fill=TRUE)

## =========================================
## 3. 분석 모델링 (Lag 0-7, SGG별 루프)
## =========================================
# 🎯 Lag를 7일까지 확장
lags <- 0:7
res_lag_list <- list()

if (!exists("target_list") || nrow(target_list) == 0) {
  all_sgg_codes <- unique(expo_all$SGG_CD)
} else {
  all_sgg_codes <- target_list$SGG_CD 
}

for (lg in lags) {
  cat("\nRunning Lag:", lg, "\n")
  
  # 3-1. Lag 변수 생성
  new_all <- copy(expo_all)[order(SGG_CD, date)]
  
  if (lg == 0) {
    new_all[, wf_lag := wf]
  } else {
    new_all[, wf_lag := shift(wf, n = lg, type = "lag"), by = SGG_CD]
  }
  
  new_all[, `:=`(wday = wday(date), month = month(date), year = year(date))]
  
  # 3-2. 데이터 합치기
  cases_only <- ha_all[respiratory > 0]
  cases_only[, id := .I] 
  
  dt_all <- merge(
    cases_only, 
    new_all, 
    by = c("SGG_CD", "wday", "month", "year"), 
    all.x = TRUE, 
    allow.cartesian = TRUE 
  )
  
  # 3-3. Case / Control 구분
  dt_all[, case := fifelse(case_date == date, 1L, 0L)]
  
  dt_all <- dt_all[!is.na(wf_lag)]
  
  # 3-4. 시군구별 루프 분석
  res_sgg <- vector("list", length(all_sgg_codes))
  
  for (i in seq_along(all_sgg_codes)) {
    sgg <- all_sgg_codes[i]
    foo <- dt_all[SGG_CD == sgg]
    
    n_case_val <- if(nrow(foo) > 0) sum(foo$case) else 0
    
    if (nrow(foo) > 0 && n_case_val > 0 && length(unique(foo$wf_lag)) >= 2) {
      
      fit <- try(
        clogit(case ~ wf_lag + strata(id), 
               data = foo, 
               weights = respiratory,
               method = "approximate"), 
        silent = TRUE
      )
      
      if (!inherits(fit, "try-error")) {
        ci <- summary(fit)$conf.int[1, c(1,3,4)]
        res_sgg[[i]] <- data.table(
          SGG_CD = sgg, lag_day = lg, 
          RR = ci[1], RR_ll = ci[2], RR_ul = ci[3], 
          n_case = n_case_val, note = "Success"
        )
      } else {
        res_sgg[[i]] <- data.table(
          SGG_CD = sgg, lag_day = lg, 
          RR = NA_real_, RR_ll = NA_real_, RR_ul = NA_real_, 
          n_case = n_case_val, note = "Model Error"
        )
      }
      
    } else {
      res_sgg[[i]] <- data.table(
        SGG_CD = sgg, lag_day = lg, 
        RR = NA_real_, RR_ll = NA_real_, RR_ul = NA_real_, 
        n_case = n_case_val, note = "Insufficient Data"
      )
    }
    
    if (i %% 50 == 0) cat(".")
  }
  
  res_lag_list[[as.character(lg)]] <- rbindlist(res_sgg)
}

# 최종 결과 합치기
res_sgg_all <- rbindlist(res_lag_list)

## =========================================
## 4. 최종 저장 (SGG 정보 포함)
## =========================================
res_all <- rbindlist(res_lag_list)
final_res <- merge(target_list, res_all, by="SGG_CD", all.x=TRUE)
final_res <- final_res[order(SIDO_NM, SGG_NM, lag_day)]

cat("Total Rows:", nrow(final_res), "\n")
print(final_res[SGG_NM %in% c("군위군", "미추홀구", "부천시", "원미구")])

fwrite(final_res, file.path(out_dir, "SGG_level_OR_2015_2024_FILLED.csv"), bom=TRUE)
cat("Saved to:", file.path(out_dir, "SGG_level_OR_2015_2024_FILLED.csv"))
