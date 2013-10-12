## imageAnalysis_function.R
## functions to accompany imageAnalysis.R
## thus far poorly commented


library(tiff)
# library(rtiff)
library(EBImage)
# library(ggplot2)

## take the projection of an arbitrary object on a channel at a given time
## and return a vector corresponding to the intensity in each slice in the
## area of the object's projection. This does NOT make any normalization
## for object area, which is something that might be interesting too
get_intensity_profile = function(obj_projection , channel , time) {
    return(apply(image_array[,,channel,time,] , MARGIN = 3 ,
                 FUN = function(image_matrix , nuc_projection) {
                     sum(obj_projection*image_matrix)} ,
                 obj_projection))
## MARGIN is 3 because the matrix being passed is projected onto channel, time
}

## this can be a lot smarter... gotta think more about scope, and which variable should be where...
## Return the intensity profile at a given time for a nucleus indexed by nucleus_index
## corresponding to a row entry in kernel_table. kernel table contains the x,y
## positions of acceptable nucleii. Nucleus_projections contains the max project
## on the nucleus channel, nucleus_labels is the standard EBImage label matrix for nucleii
get_nucleus_intensity_profile =
    function(nucleus_index , kernel_table , nucleus_projections , nucleus_labels , time) {
        seed = matrix(0,nrow = 512 , ncol = 512)
        seed[kernel_table[nucleus_index,1], kernel_table[nucleus_index,2]] = 1
        nuc_projection = propagate(x = nucleus_projections , seeds = seed , mask = nucleus_labels)
        ## ^^ this is a very useful matrix. contains positional information
        ## for the indicated nucleus, lets me multiply times the intensity
        ## of other images to obtain the intensity in that vertical stack
        get_intensity_profile(nuc_projection , 1 , time)
}

get_focus_intensity_profile =
    function(focus_index , kernel_table , focus_projections , focus_labels , time , normalize = TRUE) {
        seed = matrix(0,nrow = 512 , ncol = 512)
        seed[kernel_table[focus_index,1], kernel_table[focus_index,2]] = 1
        foc_projection = propagate(x = focus_projections , seeds = seed , mask = focus_labels)
        ## ^^ this is a very useful matrix. contains positional information
        ## for the indicated nucleus, lets me multiply times the intensity
        ## of other images to obtain the intensity in that vertical stack
        profile = get_intensity_profile(foc_projection , 2 , time)
        if(normalize == TRUE) {
            profile = profile / sum(seed)
        }
        return(profile)
    }


## return a matrix where each pixel is given as the maximum intensity in each pixel stacked
## vertically on top in a single Z-stack on a single channel
max_project = function(time , channel){
    max_output = matrix(0 , ncol = length(image_array[,1,1,1,1]) , nrow = length(image_array[1,,1,1,1]))
    for(i in 1:length(image_array[,1,1,1,1])) {
        for(j in 1:length(image_array[1,,1,1,1])){
            max_output[i,j] = max(image_array[i,j,channel,time,])
        }
    }
    return(max_output)
}

## return a matrix of mean values and sd values for the area in each frame
## within w x-pixels and h y-pixels of the centroid of the focus 
get_focus_context = function(focus_index , time_point , foci , h , w) {
    address = get_focus_address(focus_index , foci)
    context = image_array[address[1]-w:address[1]+w , address[2]-h:address[2]+h , 2,time_point,]
    mean = NULL
    for(i in 1:length(context[1,1,])) {mean = rbind(mean , mean(c(context[1:length(context[,1,1]),,i] , context[,1:length(context[1,,1]),i])))}
    devs = NULL
    for(i in 1:length(context[1,1,])) {devs = rbind(devs , var(c(context[1:length(context[,1,1]),,i] , context[,1:length(context[1,,1]),i]))^2)}
    return(cbind(mean , devs))
}

## return the pixel intensity for each pixel in every slice that falls within
## w x-pixels and h y-pixels of the centroid of the focus
get_raw_focus_context = function(focus_index , time_point , foci , h , w) {
    address = get_focus_address(focus_index , foci)
    address = floor(address)
    if((address[1] - h) < 0 | (address[1]+h) > length(image_array[,1,1,1,1]) |
       (address[2] - w) < 0 | (address[2]+2) > length(image_array[1,,1,1,1]))
        { print ('ERROR: window out of bounds; focus is too close to the endge of the image?')
          return(0)}
    context = image_array[(address[1]-h):(address[1]+h) , (address[2]-w):(address[2]+w) , 2,time_point,]
    out = NULL
    for(x in 1:length(context[,1,1])){
        for(y in 1:length(context[x,,1])){
            for(z in 1:length(context[x,y,])){
                out = rbind(out , cbind(context[x,y,z],z))
            }
        }
    }
    return(cbind(pixel_value = out[,1] , z_pos = out[,2]))
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
