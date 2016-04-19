# Logs dirs.
dir_logs = paste("logs", "/", sep = "")
dir_def  = paste(dir_logs, "default", "/", sep = "")
dir_fast = paste(dir_logs, "fast", "/", sep = "")

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

# Building plots alltogether.
# Preparing main graph with limits.
plot(x = c(1), y = def_simple$Time.ns., 
     xlab = "Количество блоков", ylab = "Время (нс)",
     xlim = range(c(0, def_blocks[nrow(def_blocks), 3])), 
     ylim = range(c(12500000000, def_simple[1, 3])), 
     type = "p")
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
legend('topright', lty = 1,  bty = 'n', cex = .75,
       c("Blocks", "Lines", "Cols", "Lines and cols", "Lines w. vector."), 
       col=c('red', 'green', 'blue', 'orange', 'purple'))

# Building only default parallel blocks.
# Preparing main graph with limits.
yMax = def_par_ml[1, 4]
if (def_par_mc[1, 4] > yMax)
  yMax = def_par_mc[1, 4]
if (def_par_mlc[1, 4] > yMax)
  yMax = def_par_mlc[1, 4]
if (def_par_mlv[1, 4] > yMax)
  yMax = def_par_mlv[1, 4]
plot(x = def_par_ml$Blocksize.m., y = def_par_ml$Time.ns., 
     xlab = "Количество блоков", ylab = "Время (нс)",
     xlim = range(c(0, def_par_ml[nrow(def_blocks), 3])), 
     ylim = range(c(29000000000, yMax)), 
     type = "l", col = "green")
# Adding all default evals.
lines(x = def_par_mc$Blocksize.m.,  y = def_par_mc$Time.ns.,  col = "blue")
lines(x = def_par_mlc$Blocksize.m., y = def_par_mlc$Time.ns., col = "orange")
lines(x = def_par_mlv$Blocksize.m., y = def_par_mlv$Time.ns., col = "red")
legend('topright', lty = 1,  bty = 'n', cex = .75,
       c("Lines", "Cols", "Lines and cols", "Lines w. vector."), 
       col=c('green', 'blue', 'orange', 'red'))

# Building only fast parallel blocks.
# Preparing main graph with limits.
yMax = fast_par_ml[1, 4]
if (fast_par_mc[1, 4] > yMax)
  yMax = fast_par_mc[1, 4]
if (fast_par_mlc[1, 4] > yMax)
  yMax = fast_par_mlc[1, 4]
if (fast_par_mlv[1, 4] > yMax)
  yMax = fast_par_mlv[1, 4]
plot(x = fast_par_ml$Blocksize.m., y = fast_par_ml$Time.ns., 
     xlab = "Количество блоков", ylab = "Время (нс)",
     xlim = range(c(0, fast_par_ml[nrow(fast_blocks), 3])), 
     ylim = range(c(9500000000, yMax)), 
     type = "l", col = "green")
# Adding all fastault evals.
lines(x = fast_par_mc$Blocksize.m.,  y = fast_par_mc$Time.ns.,  col = "blue")
lines(x = fast_par_mlc$Blocksize.m., y = fast_par_mlc$Time.ns., col = "orange")
lines(x = fast_par_mlv$Blocksize.m., y = fast_par_mlv$Time.ns., col = "red")
legend('topright', lty = 1,  bty = 'n', cex = .75,
       c("Lines", "Cols", "Lines and cols", "Lines w. vector."), 
       col=c('green', 'blue', 'orange', 'red'))
