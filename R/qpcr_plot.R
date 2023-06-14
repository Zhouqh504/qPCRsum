#' @title qpcr_plot
#' @description Visualize the qPCR result through bar plot
#'     for single gene between Control and single/different treatments.
#' @param dd.ct the output of function of Calculate_dd_ct.
#' @param gene the qPCR result of gene you want to show.
#' @param t_test_results the output of function of perform_t_tests
#' @param bar_width modify the width of the bars
#' @param errorbar_bar modify the width of the errorbar
#' @param bar_position modify the spacing between the bars
#' @param bar_colors modify the colors of the bars
#' @param p_val_labelsize modify the size of p-value
#' @param axis_text_x_size modify the size of the x-axis text
#' @param axis_text_y_size modify the size of the y-axis text
#' @param axis_title_y_size modify the size of the y-axis title
#' @param plot_width modify the width of plot
#' @param plot.height modify the height of plot
#'
#' @return Bar plot of qPCR result for single gene in different treatments.
#' @export
#' @import dplyr
#' @import ggplot2
#' @examples
#' dd.ct <- data.frame(
#'   biological_replicate = c(1, 1, 2, 2, 3, 3, 4, 4),
#'   treatment = c("Control", "Treated1", "Control", "Treated1", "Control", "Treated1", "Control", "Treated1"),
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
#' t_test_results <- data.frame(
#'   column_name = c("gene01_2_ddct", "gene02_2_ddct"),
#'   group1 = c("Control", "Control"),
#'   group2 = c("Treated1", "Treated1"),
#'   p_value = c(0.255771789, 0.004903698),
#'   significance = c("ns", "**")
#' )
#' gene <- "gene01"
#' control <- "Control"
#' qpcr_plot(
#'   dd.ct = dd.ct, gene = gene, control = control, order = NULL, t_test_results = t_test_results,
#'   y_lab = NULL, bar_size = 1, bar_width = 0.8, errorbar_width = 0.2, errorbar_size = 1, bar_position = 0.8, bar_colors = NULL, p_val_labelsize = 6,
#'   axis_text_x_size = 14, axis_text_y_size = 14, axis_title_y_size = 18,
#'   y_max = NULL, panel.border.size = 1, axis.ticks.size = 1, rotation_angle = 0, axis.text.x.dis = 0.5, axis.ticks.length = 0.15
#' )
qpcr_plot <- function(dd.ct = dd.ct, gene = gene, control = control, order = NULL, t_test_results = NULL, y_lab = NULL, bar_size = 1, bar_width = 0.8, errorbar_width = 0.2, errorbar_size = 1, bar_position = 0.8, bar_colors = NULL, p_val_labelsize = 6,
                      axis_text_x_size = 14, axis_text_y_size = 14, axis_title_y_size = 18,
                      y_max = NULL, panel.border.size = 1, axis.ticks.size = 1, rotation_angle = 0, axis.text.x.dis = 0.5, axis.ticks.length = 0.15) {
  # 获取不重复的组别
  treatments <- unique(dd.ct$treatment)


  # 获取指定范围内的列名

  selected_columns <- colnames(dd.ct)[grepl("ddct", colnames(dd.ct))]

  # 定义循环的初始为0
  gene_num <- 0

  dd.cts <- NULL # 在循环之前初始化数据框

  for (column_name in selected_columns) {
    gene_num <- gene_num + 1
    column_name <- paste0(column_name, "") # 使用带有"_2_ddct"后缀的更新列名
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
  }



  if (!is.null(t_test_results)) {
    p_val <- subset(t_test_results, column_name == paste0(gene, "_2_ddct"))
    # 筛选出 treatment 为 "treated" 的数据
    p_val <- subset(p_val, group1 == control)
  }
  mean_column_name <- paste0("mean_", gene, "_2_ddct")
  sd_column_name <- paste0("sd_", gene, "_2_ddct")
  # print(mean_column_name)
  # print(sd_column_name)

  # 根据指定的顺序重新排序 x 轴变量
  if (!is.null(order)) {
    dd.cts$treatment <- factor(dd.cts$treatment, levels = order)
  }


  # 应用自定义主题到ggplot对象#修改柱形的宽度 修改柱形的间距

  # 加柱子颜色
  if (!is.null(bar_colors)) {
    p <- ggplot(dd.cts, aes(x = treatment, y = !!sym(mean_column_name), fill = treatment, width = 0.5)) +
      scale_fill_manual(values = bar_colors) +
      theme(legend.position = "none") +
      theme_bw() +
      geom_bar(stat = "identity", position = position_dodge(width = bar_position), colour = "black", width = bar_width, size = bar_size)
  }
  if (is.null(bar_colors)) {
    p <- ggplot(dd.cts, aes(x = treatment, y = !!sym(mean_column_name), width = 0.5)) +
      theme(legend.position = "none") +
      theme_bw() +
      geom_bar(stat = "identity", position = position_dodge(width = bar_position), colour = "black", width = bar_width, size = bar_size)
  }

  # 修改柱形的宽度 修改柱形的间距

  # 加误差棒 是否可以封装成一个函数
  p <- p + geom_errorbar(aes(ymin = !!sym(mean_column_name) - !!sym(sd_column_name), ymax = !!sym(mean_column_name) + !!sym(sd_column_name)), width = errorbar_width, color = "black", size = errorbar_size, position = position_dodge(width = bar_position))


  # 加横纵坐标标签
  if (is.null(y_lab)) {
    p <- p + xlab("") + ylab(expression(paste("2(-", Delta, Delta, "Ct)")))
  } else {
    p <- p + xlab("") + ylab(y_lab)
  }

  # theme
  p <- p + theme_bw() + theme(legend.position = "none")


  p <- p + theme(
    axis.text.x = element_text(face = "bold", size = axis_text_x_size, colour = "black"),
    axis.text.y = element_text(face = "bold", size = axis_text_y_size, colour = "black"),
    axis.title.y = element_text(size = axis_title_y_size, face = "bold", colour = "black")
  )

  if (!is.null(t_test_results)) {
    # 在柱状图上仅对 treated_results 数据进行标注
    p <- p + annotate("text", x = p_val$group2, y = max(dd.cts[[mean_column_name]]) + max(dd.cts[[sd_column_name]]) + 0.5, label = p_val$significance, size = p_val_labelsize)
  }

  # print(dd.cts)

  # 根据传入的参数设置纵坐标轴范围
  if (!is.null(y_max)) {
    p <- p + ylim(0, y_max)
  }

  # 设置刻度的粗细和长度
  if (axis.ticks.length > 0) {
    axis.ticks.length <- unit(axis.ticks.length, "cm")
    p <- p + theme(axis.ticks = element_line(size = axis.ticks.size, colour = "black"), axis.ticks.length = unit(axis.ticks.length, "cm")) # 设置刻度线的粗细和颜色
  } else {
    # 如果长度小于等于零，给出警告或采取适当的处理方法
    warning("Length value must be greater than zero.")
    axis.ticks.length <- unit(0.15, "cm") # 或者使用其他默认值
  }


  # 设置图形的边框线条颜色和粗细
  p <- p + theme(panel.border = element_rect(color = "black", size = 1))

  # 设置横坐标标签的倾斜

  p <- p + theme(axis.text.x = element_text(angle = rotation_angle, vjust = axis.text.x.dis))

  # 设置面板边框大小为1
  p <- p + theme(panel.border = element_rect(size = panel.border.size))




  # 添加标题和副标题，具体参数设置了解,用户设置的基因
  p <- p +
    ggtitle(gene) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, colour = "black", face = "bold")) +
    theme(panel.grid = element_blank())

  # 图尺寸
  # 图表尺寸设置
  # if (!is.null(plot_width) && !is.null(plot_height)) {
  # p <- p +theme(plot.width = plot_width, plot.height = plot_height)
  # }

  # 返回 ggplot 图形对象
  return(p)
}
