library(dplyr)
library(readr)
library(meta)
library(stringr)
library(ggplot2) # 버블 플롯을 위한 라이브러리 추가 (bubble 함수가 ggplot2에 의존할 수 있음)

# 1. 데이터 불러오기 및 전처리
# -----------------------------------------------------------------------------
path_or <- "C:/Users/USER/Desktop/thesis/output3/SGG_level_RR_2015_2024_wildfire_affected_cc_hourly.csv"
path_comm <- "C:/Users/USER/Desktop/thesis/community/community.csv"
# df_or 파일명: SGG_level_OR_2015_2024_FILLED.csv 로 확인되었습니다.

df_or <- read_csv(path_or, show_col_types = FALSE)
df_comm <- read_csv(path_comm, show_col_types = FALSE)

vars_to_clean <- c("edu", "work", "green", "old", "poor", "wealth")

# Community 데이터 전처리
df_comm_clean <- df_comm %>%
  mutate(across(all_of(vars_to_clean), 
                ~ as.numeric(str_replace_all(., "%", "")), 
                .names = "{.col}")) %>%
  mutate(SGG_CD = as.character(SGG_CD))

# Health (OR) 데이터 전처리 & SE 계산
df_or_clean <- df_or %>%
  # ✅ 수정된 부분: lag == "same_day" 대신 lag_day == 1 사용
  filter(lag_hour == 1) %>% 
  mutate(
    logRR = log(RR),
    # 신뢰구간을 이용해 표준오차(seTE) 계산
    seTE = (log(RR_ul) - log(RR_ll)) / 3.92,
    SGG_CD = as.character(SGG_CD)
  ) %>%
  filter(!is.na(logRR) & !is.na(seTE))

# 병합
merged_data <- inner_join(df_comm_clean, df_or_clean, by = "SGG_CD")

# -----------------------------------------------------------------------------
# 2. [중요] 극단적인 분산값 확인 및 제거 (Filtering)
# -----------------------------------------------------------------------------

# 분산(Variance) 확인: seTE의 제곱
merged_data$varTE <- merged_data$seTE^2

cat("=== 필터링 전 데이터 요약 ===\n")
summary(merged_data$seTE)

# [해결책] 표준오차(seTE)가 너무 큰 값 제거
threshold_se <- 2.0  # 표준오차가 2.0 이상인 데이터는 제외 (매우 불안정한 데이터)

filtered_data <- merged_data %>%
  filter(seTE < threshold_se) %>%
  filter(seTE > 0.0001) # 혹시 모를 0 분산(너무 작은 값)도 제외

cat("\n=== 필터링 후 데이터 개수 ===\n")
cat("원본 개수:", nrow(merged_data), "-> 필터링 후:", nrow(filtered_data), "\n")
cat("제외된 시군구 수:", nrow(merged_data) - nrow(filtered_data), "\n")


# -----------------------------------------------------------------------------
# 3. 메타 회귀분석 재수행 (필터링된 데이터 사용)
# -----------------------------------------------------------------------------
results_list <- list()

# 메타분석에 사용할 변수 목록
for (var in vars_to_clean) {
  temp_data <- filtered_data %>% filter(!is.na(!!sym(var)))
  
  if(nrow(temp_data) > 3) { # 최소한의 데이터 포인트 확인 (자유도 확보)
    # method.tau를 "REML"로 변경 (불균형 데이터에서 DL보다 수렴이 잘 될 수 있음)
    m_gen <- metagen(TE = logRR,
                     seTE = seTE,
                     studlab = SGG_CD,
                     data = temp_data,
                     sm = "OR",
                     method.tau = "REML") # DerSimonian-Laird 대신 REML 사용
    
    # 메타 회귀분석
    m_reg <- metareg(m_gen, as.formula(paste("~", var)))
    
    # 결과 추출 (회귀계수, 표준오차, P-value는 두 번째 값(covariate)을 사용)
    res_row <- data.frame(
      Variable = var,
      Beta = m_reg$beta[2],
      SE = m_reg$se[2],
      P_value = m_reg$pval[2],
      Tau2 = m_reg$tau2,
      I2 = max(0, m_reg$I2)
    )
    results_list[[var]] <- res_row
  } else {
    cat(sprintf("⚠ %s 변수는 NA 제거 후 데이터가 부족하여 메타 회귀분석을 건너뜁니다. (N=%d)\n", var, nrow(temp_data)))
  }
}

final_results <- do.call(rbind, results_list)

# 결과 저장
output_path <- "C:/Users/USER/Desktop/thesis/output3/meta_regression_results.csv"
write.csv(final_results, output_path, row.names = FALSE)
print(final_results)
