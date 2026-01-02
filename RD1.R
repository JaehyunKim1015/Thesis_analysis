library(data.table)
library(lubridate)
library(survival) # survival은 RD 분석에 필요 없지만 템플릿을 위해 유지

## =========================================
## 0. 경로 설정
## =========================================
base_dir   <- "C:/Users/USER/Desktop/thesis"
# 🎯 SGG 매핑 파일 경로
map_dir    <- file.path(base_dir, "sgg_cd")
expo_dir   <- file.path(base_dir, "wf_bi")
# Health 데이터 경로 수정 (daily_count 파일명에 맞춤)
nedis_dir  <- file.path(base_dir, "NEDIS_count_daily")
pop_dir    <- file.path(base_dir, "pop_with_code")
sgg_dir    <- file.path(base_dir, "SGG")
out_rd_dir <- file.path(base_dir, "output1", "sgg_rd")

if (!dir.exists(out_rd_dir)) dir.create(out_rd_dir, recursive = TRUE)

years <- 2015:2024 # 🎯 분석 기간 2024년까지 확장

# 매핑 적용 함수 (SGG_CD 타입 통일 및 매핑)
convert_sgg <- function(dt, map_dt) {
  if ("SGG_CD" %in% names(dt)) {
    dt[, old_code := as.character(SGG_CD)]
    dt[map_dt, on = .(old_code), SGG_CD := i.new_code]
    dt[, old_code := NULL] 
    dt <- dt[!is.na(SGG_CD)]
  }
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

# 🎯 FIX: 모든 코드를 문자열로 통일
df_2022[, code := as.character(code)]
df_2023[, code := as.character(code)]
df_2024[, code := as.character(code)]

remove_space <- function(x) gsub("\\s+", "", x)
df_2022[, `:=`(sido = remove_space(sido), sigungu = remove_space(sigungu))]
df_2023[, `:=`(sido = remove_space(sido), sigungu = remove_space(sigungu))]
df_2024[, `:=`(sido = remove_space(sido), sigungu = remove_space(sigungu))]

df_hist <- rbindlist(list(
  df_2022[, .(old_code = code, sido, sigungu)],
  df_2023[, .(old_code = code, sido, sigungu)]
), use.names = TRUE)
df_hist <- unique(df_hist)

target_list <- unique(df_2024[, .(SGG_CD = code, SIDO_NM = sido, SGG_NM = sigungu)])
target_list[, SGG_CD := as.character(SGG_CD)]

base_map <- merge(
  df_hist, 
  df_2024[, .(new_code = code, sido, sigungu)], 
  by = c("sido", "sigungu"), all.x = TRUE
)

# 예외 처리
base_map[old_code == "47720" & is.na(new_code), new_code := "27720"] # 군위군

michuhol_code_2024 <- df_2024[sigungu %like% "미추홀", code][1]
if(!is.na(michuhol_code_2024)) base_map[old_code == "28170" & is.na(new_code), new_code := michuhol_code_2024] # 미추홀구

bucheon_representative_2024 <- df_2024[sido == "경기" & sigungu %like% "부천시", code][1]
bucheon_old_codes <- unique(c(
  df_hist[sigungu == "부천시", old_code], 
  df_2023[sigungu %in% c("원미구", "소사구", "오정구"), code]
))
bucheon_old_codes <- bucheon_old_codes[!bucheon_old_codes %in% df_2024[sigungu %like% "부천", code]] 

if(!is.na(bucheon_representative_2024)) {
  base_map[old_code %in% bucheon_old_codes & is.na(new_code), new_code := bucheon_representative_2024]
}

# 최종 매핑 테이블 생성
final_map_list <- list()
final_map_list[[1]] <- base_map[!is.na(old_code) & !is.na(new_code), .(old_code, new_code)]
existing_mapped_olds <- final_map_list[[1]]$old_code
self_map <- target_list[!SGG_CD %in% existing_mapped_olds, .(old_code = SGG_CD, new_code = SGG_CD)]
final_map_list[[2]] <- self_map

sgg_map <- rbindlist(final_map_list)
sgg_map <- unique(sgg_map)

cat("SGG Mapping table generated. Total unique mappings:", nrow(sgg_map), "\n")


## ==============================================================================
## 1. month.control 함수 (Schwartz et al. 원본 로직 유지)
## ==============================================================================
month.control <- function(exposed, control, nbuffer = 30) {
  if (length(exposed) > 0) {
    out <- lapply(exposed, function(e_day) {
      e_buffer <- e_day + (-nbuffer:nbuffer)
      c_pool <- control[control %in% e_buffer]
      return(c(e_day, c_pool))
    })
    names(out) <- exposed
  } else {
    out <- numeric()
  }
  return(out)
}

## ==============================================================================
## 2. month.wt.analysis 함수 (Schwartz et al. 원본 가중치 로직 적용)
## ==============================================================================
month.wt.analysis <- function(exposure, event, c.list, outcome.dt) {
  if (length(c.list) > 0) {
    out <- numeric()
    
    for (i in 1:length(c.list)) {
      dates <- c.list[[i]]
      
      if (length(dates) > 1) {
        e_day <- dates[1]
        
        c_pool <- data.table(
          date = dates,
          exposed = c(1, rep(0, length(dates) - 1))
        )
        
        diffs <- as.numeric(difftime(c_pool$date, e_day, units = "days"))
        dist_val <- round(abs(diffs), digits = 0)
        
        c_pool$wt <- ifelse(dist_val == 0, 0, 1 / dist_val)
        
        c_pool <- merge(c_pool, outcome.dt, by = "date", all.x = TRUE)
        
        c_pool[is.na(get(event)), (event) := 0]
        
        den <- sum(c_pool[exposed == 0]$wt, na.rm = TRUE)
        num <- sum(c_pool[exposed == 0, get(event)] * c_pool[exposed == 0]$wt, na.rm = TRUE)
        
        if (!is.na(den) && den > 0) {
          y_exposed <- c_pool[exposed == 1, get(event)][1]
          val <- y_exposed - (num / den)
          
          if (is.finite(val)) {
            out <- c(out, val)
          }
        }
      }
    }
    
    if (length(out) == 0) {
      temp <- simpleError("No non-infinite RD")
      temp$call <- NULL
    } else {
      temp <- c(mean(out), length(out))
    }
    
  } else {
    temp <- simpleError("No available exposed day")
    temp$call <- NULL
  }
  return(temp)
}


## ==============================================================================
## 3. 메인 실행 루프 (데이터 포맷 강제 변환 및 매핑 적용)
## ==============================================================================

bobb.control.n <- 3
nbuffer <- 30

for (YR in years) {
  cat("============================================\n")
  cat("▶", YR, "년 처리 시작\n")
  
  # -------------------------------------------------------
  # 1. 산불 데이터 (WF)
  # -------------------------------------------------------
  wf_file <- file.path(expo_dir, sprintf("%d_wfbi_with_sgg.csv", YR))
  if (!file.exists(wf_file)) { cat(" ❌ WF file not found.\n"); next }
  wf <- fread(wf_file)
  
  # SGG_CD 이름 변경 및 매핑 적용
  wf <- safe_rename_sgg(wf)
  wf <- convert_sgg(wf, sgg_map)
  wf[, SGG_CD := as.character(SGG_CD)]
  
  # 날짜 변환
  if (!inherits(wf$date, "Date")) wf[, date := lubridate::ymd(date)]
  
  wf <- wf[, .(SGG_CD, date, smoke_day)]
  
  # -------------------------------------------------------
  # 2. 응급실 데이터 (NEDIS)
  # -------------------------------------------------------
  # 🎯 파일명 수정
  nedis_file <- file.path(nedis_dir, sprintf("%d_daily_with_centroid.csv", YR))
  if (!file.exists(nedis_file)) { cat(" ❌ NEDIS file not found.\n"); next }
  ha <- fread(nedis_file, colClasses = c(SGG_CD="character")) # SGG_CD를 문자열로 로드
  
  # 컬럼명 정제
  if ("daily_count" %in% names(ha) && "date" %in% names(ha)) {
    setnames(ha, old = "daily_count", new = "respiratory")
    setnames(ha, old = "date", new = "date") # date 컬럼명 일치
  } else {
    cat(" ❌ Health file columns (daily_count/date) missing.\n"); next
  }
  
  # SGG_CD 정제 및 매핑 적용
  ha <- convert_sgg(ha, sgg_map)
  ha[, SGG_CD := as.character(SGG_CD)]
  # 날짜 정제
  if (!inherits(ha$date, "Date")) ha[, date := lubridate::ymd(ha$date)]
  
  # respiratory 컬럼 강제 정제 및 0 처리 유지
  ha[, respiratory := as.character(respiratory)]
  ha[, respiratory := gsub("[^0-9]", "", respiratory)]
  ha[, respiratory := as.numeric(respiratory)]
  ha[is.na(respiratory), respiratory := 0]
  
  ha <- ha[, .(SGG_CD, date, respiratory)]
  
  # -------------------------------------------------------
  # 3. 인구 데이터 (Pop)
  # -------------------------------------------------------
  pop_file <- file.path(pop_dir, sprintf("sgg_pop_%d_with_code.csv", YR))
  if (!file.exists(pop_file)) { cat(" ❌ Pop file not found.\n"); next }
  pop <- fread(pop_file)
  pop <- safe_rename_sgg(pop)
  pop <- convert_sgg(pop, sgg_map)
  pop[, SGG_CD := as.character(SGG_CD)]
  
  # pop 컬럼 정제 유지
  if ("pop" %in% names(pop)) {
    pop[, pop := as.character(pop)]
    pop[, pop := gsub("[^0-9.]", "", pop)]
    pop[, pop := as.numeric(pop)]
  } else {
    stop("Pop file is missing the 'pop' column.")
  }
  
  # -------------------------------------------------------
  # 4. Centroid 데이터
  # -------------------------------------------------------
  cent_file <- file.path(sgg_dir, sprintf("sgg_cd_%d_with_centroid.csv", YR))
  if (!file.exists(cent_file)) { cat(" ❌ Centroid file not found.\n"); next }
  cent <- fread(cent_file)
  cent <- safe_rename_sgg(cent)
  cent <- convert_sgg(cent, sgg_map)
  cent[, SGG_CD := as.character(SGG_CD)]
  
  # -------------------------------------------------------
  # [데이터 병합 및 검증]
  # -------------------------------------------------------
  dt <- merge(wf, ha, by = c("SGG_CD", "date"), all.x = TRUE)
  
  # ★ 중요: 병합 후 NA인 Outcome(respiratory)은 0으로 채움
  setnafill(dt, fill = 0, cols = "respiratory")
  
  # 매칭 결과 확인
  total_patients <- sum(dt$respiratory, na.rm = TRUE)
  cat(sprintf(">> 데이터 매칭 확인: %d년 총 환자 수 %d명\n", YR, total_patients))
  
  if (total_patients == 0) {
    cat("매칭 실패: 환자 수가 0명입니다. SGG_CD나 날짜 포맷을 다시 확인하세요.\n")
    next
  }
  
  # -------------------------------------------------------
  # [시군구별 루프 실행]
  # -------------------------------------------------------
  sgg_list <- sort(unique(dt$SGG_CD))
  event_name<- "respiratory" # 🎯 Outcome 변수명 통일
  exposure_var <- "smoke_day"
  
  # 결과 저장용 테이블
  bar <- data.table(
    SGG_CD = sgg_list,
    event = event_name,
    wf1 = NA_integer_,
    bobb_control_pool = NA_integer_,
    month_wt_rd_wf1 = NA_real_,
    month_wt_ngrp_wf1 = NA_integer_
  )
  
  # 실패 로그 저장용
  fail <- data.table(SGG_CD=character(), error=character())
  
  set.seed(824 + YR)
  
  for (i in seq_len(nrow(bar))) {
    sgg <- bar$SGG_CD[i]
    
    foo <- dt[SGG_CD == sgg]
    setorder(foo, date)
    
    # 1) Bobb et al. Exclusion: 노출일 근처(±3일)는 대조군에서 영구 제외
    buffer_days <- foo$date[foo$smoke_day == 1]
    
    if (length(buffer_days) > 0) {
      bobb_exclude <- unique(unlist(lapply(buffer_days, function(x) x + (-bobb.control.n:bobb.control.n))))
      day00_bobb <- foo$date[!(foo$date %in% bobb_exclude)]
    } else {
      day00_bobb <- foo$date
    }
    bar$bobb_control_pool[i] <- length(day00_bobb)
    
    # 2) Exposed Days 식별
    exposed_days <- foo$date[foo$smoke_day == 1]
    bar$wf1[i]<- length(exposed_days)
    
    # 3) Schwartz et al. Control Selection (±30일)
    baz2 <- month.control(exposed_days, day00_bobb, nbuffer)
    
    # 4) Weighted Analysis (Outcome 데이터만 뽑아서 전달)
    # 🎯 Outcome 컬럼명 통일
    outcome_subset <- foo[, .(date, respiratory)] 
    
    temp4 <- month.wt.analysis(
      exposure= exposure_var,
      event = event_name,
      c.list= baz2,
      outcome.dt = outcome_subset
    )
    
    # 5) 결과 저장
    if (!inherits(temp4, "condition")) {
      bar$month_wt_rd_wf1[i] <- temp4[1]
      bar$month_wt_ngrp_wf1[i] <- temp4[2]
    } else {
      fail <- rbind(fail, data.table(SGG_CD=as.character(sgg), error=as.character(temp4)), fill=TRUE)
      bar$month_wt_ngrp_wf1[i] <- length(baz2)
    }
    
    if (i %% 100 == 0) cat(" -", YR, "년", i, "개 시군구 처리 완료\n")
  }
  
  # -------------------------------------------------------
  # [결과 정리 및 저장]
  # -------------------------------------------------------
  # 인구 및 중심점 정보 결합
  bar <- merge(bar, pop[, .(SGG_CD, pop)], by="SGG_CD", all.x=TRUE)
  bar <- merge(bar, cent, by="SGG_CD", all.x=TRUE)
  
  # Standardized RD (인구 10만명 당) 계산
  bar[, month_wt_srd := (as.numeric(month_wt_rd_wf1) / as.numeric(pop)) * 100000]
  
  out_bar_file <- file.path(out_rd_dir, sprintf("%d_sgg_rd_monthwt.csv", YR))
  out_fail_file <- file.path(out_rd_dir, sprintf("%d_sgg_rd_monthwt_fail.csv", YR))
  
  fwrite(bar, out_bar_file, bom = TRUE)
  
  if(nrow(fail) > 0) {
    fwrite(fail, out_fail_file, bom = TRUE)
  }
  
  cat("▶", YR, "년 저장 완료:", out_bar_file, "\n")
}
