#' @title ddct_t_tests
#' @description t test is used to perform Statistical Analysis of −ΔΔCt
#'    between control and each treatment.
#' @param data A table input with special sequence columns
#'    (from left 2 right as follow:sample,treatment,biological_replicate,
#'    refer gene name,your gene name1,your gene name2,..),
#'    each column split with tab and containing header.
#' @param dd.ct the output of function of Calculate_dd_ct.
#' @param control the control condition name
#'
#' @return t_test_results
#' @export
#' @import dplyr
#' @examples
#' # Create an example data frame
#' dd.ct <- data.frame(
#'   biological_replicate = c(1, 1, 2, 2, 3, 3, 4, 4),
#'   treatment = c("Control", "Treated", "Control", "Treated", "Control", "Treated", "Control", "Treated"),
#'   gene01 = c(3.50, 2.68, 4.17, 3.11, 3.06, 3.71, 4.73, 3.67),
#'   gene02 = c(9.45, 6.85, 8.17, 5.76, 8.50, 6.80, 11.2, 7.84),
#'   treatment_control_vals = rep("Control", 8),
#'   gene01_control_vals = rep(3.50, 8),
#'   gene02_control_vals = rep(9.45, 8),
#'   gene01_diff = c(0, -0.82, 0, -1.06, 0, 0.65, 0, -1.06),
#'   gene01_2_ddct = c(1, 0.66, 1, 0.68, 1, 1.38, 1, 0.68),
#'   gene02_diff = c(0, -2.60, 0, -2.41, 0, -1.65, 0, -3.36),
#'   gene02_2_ddct = c(1, 0.73, 1, 0.70, 1, 0.81, 1, 0.70)
#' )
#'
#' # Print the example data frame
#' print(dd.ct)
#' # example
#' data <- data.frame(
#'   sample = rep(rep(1:4, each = 3), 2),
#'   treatment = rep(c("Control", "Treated"), each = 12),
#'   your_refer_gene = c(
#'     19.68, 19.69, 19.80, 19.95, 19.93, 19.97, 19.93, 20.02, 20.27,
#'     19.93, 19.88, 19.90, 20.61, 19.98, 20.57, 19.68, 19.95,
#'     19.85, 20.27, 20.08, 20.07, 20.10, 20.07, 20.10
#'   ),
#'   gene01 = c(
#'     23.22, 23.34, 23.12, 24.06, 24.15, 24.15, 23.18, 23.13,
#'     23.10, 24.78, 24.45, 24.67, 23.11, 22.99, 23.10, 22.77,
#'     22.99, 23.06, 23.73, 24.01, 23.80, 23.73, 23.83, 23.73
#'   ),
#'   gene02 = c(
#'     29.08, 29.04, 29.39, 28.23, 28.01, 28.12, 28.79, 28.43,
#'     28.49, 31.37, 30.74, 31.09, 27.11, 27.24, 27.37, 25.52,
#'     25.72, 25.52, 27.43, 26.73, 26.65, 27.96, 27.84, 27.98
#'   )
#' )
#' control <- "Control"
#' # t-test
#' t_test_results <- ddct_t_tests(data, dd.ct, control = "Control")
ddct_t_tests <- function(data, dd.ct, control) {
  # 获取不重复的组别
  treatments <- unique(dd.ct$treatment)

  # 创建一个空的列表，用于存储显著性检验结果
  t_test_results <- data.frame(column_name = character(), group1 = character(), group2 = character(), p_value = numeric())

  # 指定范围的列索引
  start_index <- 5
  end_index <- dim(data)[2]

  # 获取指定范围内的列名
  selected_columns <- names(data)[start_index:end_index]

  # 定义循环的初始为0
  gene_num <- 0

  dd.cts <- NULL # 在循环之前初始化数据框

  for (column_name in selected_columns) {
    gene_num <- gene_num + 1
    column_name <- paste0(column_name, "_2_ddct") # 使用带有"_2_ddct"后缀的更新列名
    # print(column_name)

    if (gene_num == 1) {
      dd.cts <- dd.ct %>%
        group_by(treatment) %>%
        summarize(
          !!paste("mean_", column_name, sep = "") := mean(get(column_name)),
          !!paste("sd_", column_name, sep = "") := sd(get(column_name))
        )
    } else {
      ddcts_value <- dd.ct %>%
        group_by(treatment) %>%
        summarize(
          !!paste("mean_", column_name, sep = "") := mean(get(column_name)),
          !!paste("sd_", column_name, sep = "") := sd(get(column_name))
        )

      dd.cts <- merge(dd.cts, ddcts_value, by = "treatment")
    }

    # 循环进行两两独立样本 t 检验


    treatment2 <- subset(treatments, treatments != control)
    group1 <- control
    for (i in 1:(length(treatment2))) {
      group2 <- treatment2[i]

      if (group1 == control) {
        # 执行独立样本 t 检验
        subset_data <- subset(dd.ct, dd.ct$treatment %in% c(group2))
        t_result <- t.test(log(subset_data[[column_name]], base = 2), mu = 0)
      }


      # 将结果存储在列表中，使用组合的名称作为列表的名称
      result_name <- paste(group1, group2, sep = "_vs_")
      result_name <- paste(column_name, result_name, sep = "_")
      # print(result_name)

      # 创建临时数据框来存储当前循环的结果
      temp_df <- data.frame(column_name = column_name, group1 = group1, group2 = group2, p_value = t_result[[3]])

      # 判断是否有一个组等于 "control"，如果有则添加到结果数据框中
      if (group1 == control) {
        t_test_results <- rbind(t_test_results, temp_df)
      }
    }
  }

  # 添加星号表示显著性水平
  t_test_results$significance <- ifelse(t_test_results$p_value < 0.001, "***",
    ifelse(t_test_results$p_value < 0.01, "**",
      ifelse(t_test_results$p_value < 0.05, "*", "ns")
    )
  )


  return(t_test_results)
}
