## 
## imageAnalysis.R
## Gavin Schlissel
## Levine Lab UC Berkeley
## Analyze 3-D image data for nascent transcript project
##

library(tiff)
# library(rtiff)
library(EBImage)

setwd('~/Documents/school/Grad/labWork/LevineLab/image_data/SnaESnaP_R_playpen/')

all_images = readTIFF('SnaESnaPB_RGB.tif' , native = F , all = T , info = T , indexed = T)

struct_image = list()
for(i in 1:15) {
  stack = seq(i,length(all_images) , 15)
  struct_image[[i]] = all_images[stack]
}


image_unlisted = unlist(struct_image)
image_array = array(image_unlisted , dim = c(512,512,3,34,15)) ##SLOW
# save(image_array , file = 'image_matrix.rdata')
load('image_matrix.rdata')




imageStack = get_stack(struct_image , 30)
display(imageStack)

my_image = struct_image[[8]][[30]]
my_image = image_array[,,,30,8]
foci = label_foci(my_image[,,2])
foci_features = computeFeatures(foci , my_image[,,2])

plot(get_focus_intensity(74 , 29 , foci_features) , col = 1 , cex = 0.3 , type ='l' , ylim = c(0,1))
lines(get_focus_intensity(74 , 30 , foci_features) , col = 2 , cex = 0.3)
lines(get_focus_intensity(74 , 31 , foci_features) , col = 3 , cex = 0.3)
lines(get_focus_intensity(74 , 32 , foci_features) , col = 3 , cex = 0.3)
lines(get_focus_intensity(74 , 33 , foci_features) , col = 3 , cex = 0.3)

tmp = get_focus_context(74 , 30 , foci_features, 3 , 3)

get_focus_context = function(focus_index , time_point , foci , h , w) {
    address = get_focus_address(focus_index , foci)
    context = image_array[address[1]-h:address[1]+h , address[2]-w:address[2]+w , 2,time_point,]
    mean = NULL
    for(i in 1:length(context[1,1,])) {mean = rbind(mean , mean(c(context[1:length(context[,1,1]),,i] , context[,1:length(context[1,,1]),i])))}
    devs = NULL
    for(i in 1:length(context[1,1,])) {devs = rbind(devs , var(c(context[1:length(context[,1,1]),,i] , context[,1:length(context[1,,1]),i]))^2)}
    return(cbind(mean , devs))
}

## requires the image matrix as as global variable called image_array
get_focus_intensity = function(focus_index , time_point ,  foci){
    address = get_focus_address(focus_index , foci)
    return(image_array[address[1] , address[2] , 2 , time_point ,])
}


get_focus_address = function(focus_index , foci) {
   foci[focus_index , c('x.0.m.cx' , 'x.0.m.cy')]
}

get_stack = function(imageStruct , time){
  images = lapply(1:15 , function(x) {process_image(imageStruct[[x]][[time]])})
  tmp = combine(images)
  return(tmp)
}

gaussian_filter = function(image) {
  y = filter2(image , makeBrush(7 , shape = 'gaussian' , step = F, sigma = 1))
  x = filter2(image , makeBrush(7 , shape = 'gaussian' , step = F , sigma = 1.5))
  return(y - x)
}

process_image = function(my_image) {
  img = rgbImage(green = my_image[,,2] , red = my_image[,,1])
  ## identify nucleii
  ## basic thresholding and shape specification
  nuc = my_image[,,1]
  nuc = label_nuc(nuc)
  
  ## identify nucleii by propogating dense centers of nucleii
  ## think about a way to replace the 0.2 constant!!!
  ctmask = opening(my_image[,,1] > 0.2 , makeBrush(5, shape = 'disc'))
  nuc = propagate(my_image[,,1] , seeds = nuc , mask = ctmask)
  
  ## reassemble image and plot with border around the nucleii
  res = paintObjects(nuc , img , col = c('red') , opac=c(1,0.5))
  
  ## identify foci
  foci = my_image[,,2]
  foci = label_foci(foci)
  
  ## paint onto image of nucleii
  return(paintObjects(foci , res , col = c('purple') , opac=c(1,0.05)))
}


label_foci = function(fociImage) {
  fociImage = gaussian_filter(fociImage)
  fociImage = thresh(fociImage , h = 20 , w = 20 , offset = 0.0425)
  fociImage = opening(fociImage, makeBrush(3 , shape = 'gaussian'))
  fociImage = fillHull(fociImage)
  return(bwlabel(fociImage))
}
label_nuc = function(nucImage) {
  nucImage = thresh(nucImage , h = 20 , w = 20 , offset = 0.05)
  nucImage = opening(nucImage, makeBrush(3 , shape = 'gaussian'))
  nucImage = fillHull(nucImage)
  return(bwlabel(nucImage))
}
