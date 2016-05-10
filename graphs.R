# Logs dirs.
dir_logs = paste("logs", "/", sep = "")
dir_def  = paste(dir_logs, "default", "/", sep = "")
dir_fast = paste(dir_logs, "fast", "/", sep = "")

# Combined output data file.
output_file = paste(dir_logs, "combined.csv", sep = "")
  
# Reading data.
# Default.
def_simple  = read.csv(paste(dir_def, "simple.csv", sep = ""))
def_blocks  = read.csv(paste(dir_def, "blocks.csv", sep = ""))
def_par_ml  = read.csv(paste(dir_def, "blocks_p_main_lines.csv", sep = ""))
def_par_mc  = read.csv(paste(dir_def, "blocks_p_main_cols.csv", sep = ""))
def_par_mlc = read.csv(paste(dir_def, "blocks_p_main_lines_cols.csv", sep = ""))
def_par_mlv = read.csv(paste(dir_def, "blocks_p_main_lines_vector.csv", sep = ""))
# Fast.
fast_simple  = read.csv(paste(dir_fast, "simple.csv", sep = ""))
fast_blocks  = read.csv(paste(dir_fast, "blocks.csv", sep = ""))
fast_par_ml  = read.csv(paste(dir_fast, "blocks_p_main_lines.csv", sep = ""))
fast_par_mc  = read.csv(paste(dir_fast, "blocks_p_main_cols.csv", sep = ""))
fast_par_mlc = read.csv(paste(dir_fast, "blocks_p_main_lines_cols.csv", sep = ""))
fast_par_mlv = read.csv(paste(dir_fast, "blocks_p_main_lines_vector.csv", sep = ""))

# Converting to seconds.
ns = 1000000000
def_simple[1,3] = def_simple[1,3] / ns
fast_simple[1,3] = fast_simple[1,3] / ns
for (i in 1 : nrow(def_blocks))
{
  def_blocks[i,4]  = def_blocks[i,4] / ns
  def_par_ml[i,4]  = def_par_ml[i,4] / ns
  def_par_mc[i,4]  = def_par_mc[i,4] / ns
  def_par_mlc[i,4] = def_par_mlc[i,4] / ns
  def_par_mlv[i,4] = def_par_mlv[i,4] / ns
  fast_blocks[i,4]  = fast_blocks[i,4] / ns
  fast_par_ml[i,4]  = fast_par_ml[i,4] / ns
  fast_par_mc[i,4]  = fast_par_mc[i,4] / ns
  fast_par_mlc[i,4] = fast_par_mlc[i,4] / ns
  fast_par_mlv[i,4] = fast_par_mlv[i,4] / ns
}

# Combining all data into the one table and writing results into the file.
combined_data = data.frame(matrix(NA, nrow = nrow(def_blocks), ncol = 11))
colnames(combined_data) <- c("Размер блока (m)",
                             "Блочное", "Расп. по стр.", "Расп. по столб.", "Расп. по стр. и столб.", "Расп. по стр. с вект.",
                             "Блочное -O3", "Расп. по стр. -O3", "Расп. по столб. -O3", "Расп. по стр. и столб. -O3", "Расп. по стр. с вект. -O3")
for (i in 1 : nrow(def_blocks))
{
  combined_data[i,1]  = def_blocks[i,3]
  combined_data[i,2]  = round( def_blocks[i,4], digits = 3)
  combined_data[i,3]  = round( def_par_ml[i,4], digits = 3)
  combined_data[i,4]  = round( def_par_mc[i,4], digits = 3)
  combined_data[i,5]  = round( def_par_mlc[i,4], digits = 3)
  combined_data[i,6]  = round( def_par_mlv[i,4], digits = 3)
  combined_data[i,7]  = round( fast_blocks[i,4], digits = 3)
  combined_data[i,8]  = round( fast_par_ml[i,4], digits = 3)
  combined_data[i,9]  = round( fast_par_mc[i,4], digits = 3)
  combined_data[i,10] = round( fast_par_mlc[i,4], digits = 3)
  combined_data[i,11] = round( fast_par_mlv[i,4], digits = 3)
}
write.table(combined_data, file = output_file, append = FALSE, 
            quote = FALSE, sep = ",", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE, col.names = TRUE, 
            qmethod = c("escape", "double"), fileEncoding = "")



# Building plots alltogether.
# Preparing main graph with limits.
plot(x = c(1), y = def_simple$Time.ns., 
     xlab = "Размер блока (m)", ylab = "Время вычисления произведения матриц (сек)",
     xlim = range(c(0, def_blocks[nrow(def_blocks), 3])), 
     ylim = range(c(0, 150)), 
     type = "l")
# Adding all default evals.
lines(x = def_blocks$Blocksize.m.,  y = def_blocks$Time.ns.,  col = "red")
lines(x = def_par_ml$Blocksize.m.,  y = def_par_ml$Time.ns.,  col = "green")
lines(x = def_par_mc$Blocksize.m.,  y = def_par_mc$Time.ns.,  col = "blue")
lines(x = def_par_mlc$Blocksize.m., y = def_par_mlc$Time.ns., col = "orange")
lines(x = def_par_mlv$Blocksize.m., y = def_par_mlv$Time.ns., col = "purple")
# Adding all O3 evals.
lines(x = fast_blocks$Blocksize.m.,  y = fast_blocks$Time.ns.,  col = "red")
lines(x = fast_par_ml$Blocksize.m.,  y = fast_par_ml$Time.ns.,  col = "green")
lines(x = fast_par_mc$Blocksize.m.,  y = fast_par_mc$Time.ns.,  col = "blue")
lines(x = fast_par_mlc$Blocksize.m., y = fast_par_mlc$Time.ns., col = "orange")
lines(x = fast_par_mlv$Blocksize.m., y = fast_par_mlv$Time.ns., col = "purple")
legend('topright', lty = 1,  bty = 'n', cex = .5,
       c("Блочное вычисление", "Параллельное вычисление по строкам", 
         "Параллельное вычисление по стоблцам", "Параллельное вычисление по строкам и столбцам", 
         "Параллельное вычисление по строкам с векторизацией"), 
       col=c('red', 'green', 'blue', 'orange', 'purple'))

# Building only default parallel blocks.
# Preparing main graph with limits.
plot(x = def_par_ml$Blocksize.m., y = def_par_ml$Time.ns., 
     xlab = "Размер блока (m)", ylab = "Время вычисления произведения матриц (сек)",
     xlim = range(c(0, def_par_ml[nrow(def_blocks), 3])), 
     ylim = range(c(29, 65)), 
     type = "l", col = "green")
# Adding all default evals.
lines(x = def_par_mc$Blocksize.m.,  y = def_par_mc$Time.ns.,  col = "blue")
lines(x = def_par_mlc$Blocksize.m., y = def_par_mlc$Time.ns., col = "orange")
lines(x = def_par_mlv$Blocksize.m., y = def_par_mlv$Time.ns., col = "red")
legend('topright', lty = 1,  bty = 'n', cex = .6,
       c("Параллельное вычисление по строкам", 
         "Параллельное вычисление по стоблцам", "Параллельное вычисление по строкам и столбцам", 
         "Параллельное вычисление по строкам с векторизацией"), 
       col=c('green', 'blue', 'orange', 'red'))

# Building only fast parallel blocks.
# Preparing main graph with limits.
plot(x = fast_par_ml$Blocksize.m., y = fast_par_ml$Time.ns., 
     xlab = "Размер блока (m)", ylab = "Время вычисления произведения матриц (сек)",
     xlim = range(c(0, fast_par_ml[nrow(fast_blocks), 3])), 
     ylim = range(c(9, 40)), 
     type = "l", col = "green")
# Adding all fastault evals.
lines(x = fast_par_mc$Blocksize.m.,  y = fast_par_mc$Time.ns.,  col = "blue")
lines(x = fast_par_mlc$Blocksize.m., y = fast_par_mlc$Time.ns., col = "orange")
lines(x = fast_par_mlv$Blocksize.m., y = fast_par_mlv$Time.ns., col = "red")
legend('topright', lty = 1,  bty = 'n', cex = .475,
       c("Параллельное вычисление по строкам с директивой -O3", 
         "Параллельное вычисление по стоблцам с директивой -O3", 
         "Параллельное вычисление по строкам и столбцам с директивой -O3", 
         "Параллельное вычисление по строкам с векторизацией с директивой -O3"), 
       col=c('green', 'blue', 'orange', 'red'))
