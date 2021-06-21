if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
# This code computes mae and r2 for pyseer predictions

viewop = function(bla, nm = "out") {
  print(paste(nm, ": ", round(mean(bla),3), " (", round(sd(bla),3), "), max:", max(bla), ", min:", min(bla), sep = ""))
}

fldr = 'LOSO'
pheno = "cd"
files = dir(fldr)
files = files[grep(pheno, files)]

files = files[-grep("tst.samples", files)]

pred_files = sort(files[grep("preds", files)])
tru_files = sort(files[grep("tst", files)])

res = list()

for (i in 1:length(pred_files)) {
  temp_pred = read.table(file.path(fldr, pred_files[i]), header = T)
  temp_truth = read.table(file.path(fldr, tru_files[i]), header = T)
  res[[i]] = data.frame(sample_pred = temp_pred$Sample, prediction = temp_pred$Prediction, sample_true = temp_truth$samples, truth = temp_truth$ph)
}

mae = c()
mse = c()
r2 = c()

for(i in 1:length(res)){
  print(paste("Alignment check?", all(res[[i]]$sample_pred == res[[i]]$sample_true)))
  mae = c(mae, mean(abs(res[[i]]$prediction - res[[i]]$truth)))
  mse = c(mse, mean((res[[i]]$prediction - res[[i]]$truth)^2))
  r2 = c(r2, cor(res[[i]]$prediction, res[[i]]$truth))
}
viewop(mse, "mse")
viewop(r2, "cor")
# variants used xfcv
# cd:      3454, 2725, 3667, 2766, 3258, 3383, 3093, 2952, 2912, 3159
# pen.mic: 2793, 2739, 2677, 2739, 2712, 2743, 2705, 2720, 2699, 2813
# cef.mic: 2541, 2441, 2617, 2512, 2468, 2433, 2531, 2453, 2499, 2427

# variants used loso
# cd:      3131, 3185, 3302, 3313, 3251, 3259, 3536, 3391, 3312, 3386, 3237, 3218, 3398, 3484, 3015, 3755, 3654, 3476, 3391, 3427, 3127, 3363, 3663, 3398, 3350, 3048, 3036, 3173
# pen.mic: 2675, 2639, 2678, 2665, 2736, 2706, 2731, 2745, 2666, 2654, 2689, 2650, 2745, 2659, 2573, 2874, 2604, 2715, 2568, 2407, 2713, 2720, 2626
# cef.mic: 2531, 2433, 2487, 2496, 2502, 2502, 2490, 2489, 2469, 2531, 2548, 2439, 2416, 2559, 2415, 2764, 2344, 2438, 2385, 2203, 2487, 2391, 2382