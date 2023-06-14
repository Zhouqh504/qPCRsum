#' @title calculate_2_ddct
#' @description The function is used for calculating the
#'     2(-delta delta CT) and the result is stored in a column
#'     named containing "_2_ddct".
#'
#' @param data A table input with special sequence columns
#'     (from left 2 right as follow:sample,treatment,biological_replicate,
#'     refer gene name,your gene name1,your gene name2,..),
#'     each column split with tab and containing header.
#' @param ref_gene your refer gene name
#' @param control_name the control condition name in the treatment column
#'     of data.
#'
#'
#' @return the result of 2(-delta delta CT).
#' @export
#' @import dplyr
#' @examples
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
#'
#' # Using calculate_dd_ct calculating 2(-delta delta CT)
#' dd.ct <- calculate_2_ddct(data, "your_refer_gene", "Control")
#' head(dd.ct)
#' # Values in gene01_2_ddct and gene02__2_ddct are 2(-delta delta CT)
calculate_2_ddct <- function(data, ref_gene, control_name) {
  # 去空白行
  data <- data[!apply(data, 1, function(row) all(is.na(row) | row == "")), ]

  # 列名替换
  names(data)[3] <- "biological_replicate"
  names(data)[2] <- "treatment"
  names(data)[4] <- ref_gene

  ### 作差计算
  num <- 0
  # 指定范围的列索引
  start_index <- 4
  end_index <- dim(data)[2]
  result <- NULL
  # 获取指定范围内的列名
  selected_columns <- names(data)[start_index:end_index]

  for (column_name in selected_columns) {
    num <- num + 1
    if (num == 1) {
      # 按照样本、基因和生物学重复分组，计算列的均值
      result <- data %>%
        group_by(biological_replicate, treatment) %>%
        summarise(!!column_name := mean(!!sym(column_name), na.rm = TRUE))
    } else {
      result_temp <- data %>%
        group_by(biological_replicate, treatment) %>%
        summarise(!!column_name := mean(!!sym(column_name), na.rm = TRUE))
      result <- cbind(result, result_temp[, 3])
    }
  }

  ## Delta计算
  num <- 0
  start_index <- 5
  end_index <- dim(data)[2]
  # 获取指定范围内的列名
  selected_columns <- names(data)[start_index:end_index]

  for (column_name in selected_columns) {
    num <- num + 1
    if (num == 1) {
      # 按照样本、基因和生物学重复分组，计算ct
      delta.ct <- result %>%
        group_by(biological_replicate, treatment) %>%
        summarise(!!column_name := !!sym(column_name) - !!sym(ref_gene))
    } else {
      delta.ct_temp <- result %>%
        group_by(biological_replicate, treatment) %>%
        summarise(!!column_name := !!sym(column_name) - !!sym(ref_gene))
      delta.ct <- cbind(delta.ct, delta.ct_temp[, 3])
    }
  }

  rm(delta.ct_temp)
  rm(result_temp)

  # 计算dd.ct
  control_vals <- delta.ct %>%
    filter(treatment == control_name)

  dd.ct <- delta.ct %>%
    left_join(control_vals, by = "biological_replicate", suffix = c("", "_control_vals"))

  gene_cols <- names(data)[start_index:end_index]

  for (col in gene_cols) {
    diff_col <- paste0(col, "_diff")
    ratio_col <- paste0(col, "_2_ddct")

    dd.ct <- dd.ct %>%
      group_by(biological_replicate) %>%
      mutate(
        {{ diff_col }} := ifelse(treatment == control_name, 0, get(col) - get(paste0(col, "_control_vals"))),
        {{ ratio_col }} := ifelse(treatment == control_name, 1, 2^(-1 * (get(col) - get(paste0(col, "_control_vals")))))
      )
  }

  return(dd.ct)
}
