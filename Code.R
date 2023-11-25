
# code by L.A. Courtenay
# CNRS, PACEA UMR 5199, University of Bordeaux
# last update 25/07/2023

# libraries and functions ------------

for (required_library in c("spatstat", "argosfilter", "vegan", "sf", "spData", "sp")) {
  
  # spatstat for the simulation of spatial data
  # argosfilter for some of the distance calculations
  # vegan for mantel test and jacquard distance calculations
  # sf misc library for spatial data
  # spData misc library for spatial data
  # sp misc library for spatial data
  
  if(required_library %in% rownames(installed.packages()) == TRUE) {
    
    print(paste("The package named", required_library, "has already been installed", sep = " "))
    
  } else {
    
    install.packages(required_library)
    
  }
  
}; rm(required_library)

library(spatstat)
library(sp)
library(spData)
library(vegan)
library(sf)
library(argosfilter)

# function for the calculation of parameters used to align two coordinate systems
# using the helmert transformation

helmert_params <- function(local_points, geo_points) {
  
  x_centre_local <- mean(t(local_points)[1,])
  y_centre_local <- mean(t(local_points)[2,])
  x_centre_geo <- mean(t(geo_points)[1,])
  y_centre_geo <- mean(t(geo_points)[2,])
  tx <- x_centre_geo - x_centre_local
  ty <- y_centre_geo - y_centre_local
  az_A_local <- atan2((local_points[1,1] - x_centre_local), (local_points[1,2] - y_centre_local))
  az_A_geo <- atan2((geo_points[1,1] - x_centre_geo), (geo_points[1,2] - y_centre_geo))
  az_B_local <- atan2((local_points[2,1] - x_centre_local), (local_points[2,2] - y_centre_local))
  az_B_geo <- atan2((geo_points[2,1] - x_centre_geo), (geo_points[2,2] - y_centre_geo))
  az_A <- az_A_local - az_A_geo
  az_B <- az_B_local - az_B_geo
  az <- (az_A + az_B) / 2
  
  return(list(
    x_centre_local = x_centre_local,
    y_centre_local = y_centre_local,
    azimuth = az[[1]],
    trans_x = tx,
    trans_y = ty
  ))
  
}

# function for performing the helmert transformation

helmert_transform <- function(points, helmert_params) {
  
  transformed_points <- points
  for (point in 1:nrow(points)) {
    point_centred <- c(
      (points[point,1] - helmert_params$x_centre_local),
      (points[point,2] - helmert_params$y_centre_local)
    )
    point_rotated <- c(
      (point_centred[1] * cos(helmert_params$azimuth) - point_centred[2] * sin(helmert_params$azimuth)),
      (point_centred[1] * sin(helmert_params$azimuth) + point_centred[2] * cos(helmert_params$azimuth))
    )
    point_transformed <- c(
      (point_rotated[1] + helmert_params$x_centre_local + helmert_params$trans_x),
      (point_rotated[2] + helmert_params$y_centre_local + helmert_params$trans_y)
    )
    transformed_points[point,] <- point_transformed
  }

  return(transformed_points)
  
}

# function to calculate the great circle distance

great_arc_distance <- function(coordinates) {
  
  res<-matrix(0, nrow = dim(coordinates)[1], ncol = dim(coordinates)[1])
  
  for (i in 1:dim(coordinates)[1]) {
    
    for (j in 1:dim(coordinates)[1]) {
      
      if(i != j) {
      
        res[i, j] <- argosfilter::distance(
          lat1 = coordinates[i, 2],
          lat2 = coordinates[j, 2],
          lon1 = coordinates[i, 1],
          lon2 = coordinates[j, 1]
        )
        
      }
      
    }
    
  }
  
  res[is.nan(res)] = 0
  
  return(as.dist(res))
  
}

# code to simulate random data within a predefined spatial window

simulate_random_spatial_data <- function(spatial_window, sample_size, plot = FALSE) {
  
  spatial_window <- spatstat.geom::as.owin(sf::st_cast(sf::st_multipoint(as.matrix(spatial_window)), "POLYGON"))
  
  sampled_data <- spatstat.random::rpoispp(10, win = spatial_window)
  
  coordinates <- cbind(sampled_data$x, sampled_data$y)
  coordinates <- coordinates[sample(nrow(coordinates), sample_size),]
  
  if (plot) {
    plot(spatial_window, main = "")
    points(coordinates, pch = 19)
  }
  
  return(coordinates)
  
}

# code to simulate random data within a predefined region
# spatial_window - the spatial window we are using for the overall study
# region - the region we want to simulate data within
# crs - the coordinate reference system we are using for the spData world database
# sample_size - the number of points to be simulated within this region

simulate_regional_spatial_data <- function(spatial_window, region, crs, sample_size) {
  
  # please note that values are scaled approximately to fit similar values to a long-lat coordinate system
  
  data(world)
  
  spatial_window <- sf::st_cast(sf::st_multipoint(as.matrix(spatial_window)), "POLYGON")
  
  regional_coordinates <- sf::st_as_sf(world["geom"][world$name_long == region,]) %>%
    sf::st_transform(crs = crs)
  
  simulated_coordinates <- spatstat.random::rpoispp(0.000000001, win = spatstat.geom::as.owin(regional_coordinates))
  
  # ensure that the coordinates from the spData world dataset coincide with our spatial window of 
  # mainland Europe
  
  simulated_coordinates_clean <- cbind(simulated_coordinates$x, simulated_coordinates$y)[
    spatstat.geom::inside.owin(
      simulated_coordinates$x / 100000, simulated_coordinates$y / 100000, spatstat.geom::as.owin(spatial_window)
    ),
  ]
  
  indices <- sample(nrow(simulated_coordinates_clean), sample_size, replace = FALSE)
  output_data <- simulated_coordinates_clean[indices,] / 100000
  return(output_data)
  
}

# function to simulate a matrix containing binary values describing the presence or absence of beads
# probability_vector - the probability of drawing from a bag of beads a bead of a certain type

create_attribute_matrix <- function(probability_vector,
                                    num_sites,
                                    num_bead_types,
                                    num_beads) {
  
  if (length(probability_vector) != num_bead_types) {
    stop("Probability vector has to be the same length as the number of types of beads")
  }
  
  if (sum(probability_vector) != 1) {
    stop("Probability vector must add up to 1")
  }
  
  output_matrix <- array(numeric(), dim = c(num_sites, num_bead_types))
  
  bag_of_beads <- c()
  
  for (bead in 1:length(probability_vector)) {
    
    bag_of_beads <- c(bag_of_beads, rep(bead, probability_vector[bead] * 1000))
    
  }
  
  for (site in 1:num_sites) {
    
    single_handful <- c()
    
    while(length(single_handful) != num_beads) {
      
      single_handful <- bag_of_beads[sample(1000, num_beads)]
      single_handful <- unique(single_handful)
      
    }
    
    single_handful_binary <- c()
    
    for (bead_type in 1:num_bead_types) {
      
      if (bead_type %in% single_handful) {
        single_handful_binary <- c(single_handful_binary, 1)
      } else {
        single_handful_binary <- c(single_handful_binary, 0)
      }
      
    }
    
    output_matrix[site,] <- single_handful_binary
    
  }
  
  return(output_matrix)
  
}

#

# Import Spatial Window -------------------------

# this is the coordinate reference system that was used to extract the outline of mainland europe

crs <- "+proj=laea +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

# load outline of European mainland (file without islands)

europe <- read.table(
  "Europe Coordinate File.csv", sep = ",", head = FALSE
) * 3846.125400 # file is recorded with pixel values so we scale to CRS coordinate values

# georeference and imported align the outline of europe with the spData world databases
# map of Europe

pointA <- c(632, 593) * 3846.125400
pointB <- c(799, 87) * 3846.125400
pointC <- c(1568, 627) * 3846.125400
pointA_CRS <- c(777049.180327869, 5970977.135461605)
pointB_CRS <- c(1422950.8196721314, 4025733.3908541845)
pointC_CRS <- c(4390163.9344262285, 6117730.802415876)

local_points <- rbind(pointA, pointB, pointC)
geo_points <- rbind(pointA_CRS, pointB_CRS, pointC_CRS)

helmert_parameters <- helmert_params(local_points, geo_points)
projected_europe <- helmert_transform(
  europe, helmert_parameters
) / 100000 # scale these values so they have the same approximate scale as long-lat coordinates

#

# model 1 - no spatial correlation or attribute correlation -------------

# set parameters

num_bead_types = 135
num_beads = 4
num_sites = 114
num_groups = 3
sample_size = num_sites / num_groups

# simulate spatial data

group1_coords <- simulate_random_spatial_data(projected_europe, sample_size)
group2_coords <- simulate_random_spatial_data(projected_europe, sample_size)
group3_coords <- simulate_random_spatial_data(projected_europe, sample_size)

#visualise spatial data

plot(projected_europe, main = " ", asp = 1, type = "l")
points(group1_coords[,1], group1_coords[,2], pch = 19, cex = 1.5, col = "black")
points(group2_coords[,1], group2_coords[,2], pch = 19, cex = 1.5, col = "red")
points(group3_coords[,1], group3_coords[,2], pch = 19, cex = 1.5, col = "#4d72b2")

# simulate matrices for each group

group1_beads <- create_attribute_matrix(
  rep((1 / 135), 135), # the probability of drawing any type of bead is the same
  num_sites = sample_size,
  num_bead_types = num_bead_types,
  num_beads = num_beads
)
group2_beads <- create_attribute_matrix(
  rep((1 / 135), 135), # the probability of drawing any type of bead is the same
  num_sites = sample_size,
  num_bead_types = num_bead_types,
  num_beads = num_beads
)
group3_beads <- create_attribute_matrix(
  rep((1 / 135), 135), # the probability of drawing any type of bead is the same
  num_sites = sample_size,
  num_bead_types = num_bead_types,
  num_beads = num_beads
)

atrribute_matrix <- rbind(group1_beads, group2_beads, group3_beads)
write.table(atrribute_matrix, "model1_binary_attributes.csv", row.names = FALSE, col.names = FALSE, sep = ",")
coordinate_matrix <- rbind(group1_coords, group2_coords, group3_coords)
write.table(coordinate_matrix, "model1_coordinates.csv", row.names = FALSE, col.names = FALSE, sep = ",")

# calculate disimilarity and geographical distances

jacquard_attributes <- as.matrix(vegan::vegdist(atrribute_matrix, method='jaccard'), binary = TRUE)
coordinate_distances <- great_arc_distance(coordinate_matrix)

# calculate and export mantel test results

sink("Model 1 Mantel Results.txt") # export mantel test results to matrix to txt file
vegan::mantel(jacquard_attributes, coordinate_distances, permutations = 9999)
sink()

#

# model 2 - slight overlap between groups -------------

# set parameters

num_bead_types = 135
num_beads = 4
num_sites = 114
num_groups = 3
sample_size = num_sites / num_groups

# simulate spatial data

group1_coords <- rbind(
  simulate_regional_spatial_data(projected_europe, "Spain", crs, 30),
  simulate_regional_spatial_data(projected_europe, "France", crs, 8)
)
group2_coords <- rbind(
  simulate_regional_spatial_data(projected_europe, "France", crs, 30),
  simulate_regional_spatial_data(projected_europe, "Spain", crs, 4),
  simulate_regional_spatial_data(projected_europe, "Germany", crs, 4)
)
group3_coords <- rbind(
  simulate_regional_spatial_data(projected_europe, "France", crs, 8),
  simulate_regional_spatial_data(projected_europe, "Germany", crs, 30)
)

#visualise spatial data

plot(projected_europe, main = " ", asp = 1, type = "l")
points(group1_coords[,1], group1_coords[,2], pch = 19, cex = 1.5, col = "black")
points(group2_coords[,1], group2_coords[,2], pch = 19, cex = 1.5, col = "red")
points(group3_coords[,1], group3_coords[,2], pch = 19, cex = 1.5, col = "#4d72b2")

# simulate matrices for each group

group1_beads <- create_attribute_matrix(
  c(
    rep((0.9 / 45), 45), # 90% more likely to pick up certain beads than others
    rep((0.05 / 45), 45),
    rep((0.05 / 45), 45)
  ),
  num_sites = sample_size,
  num_bead_types = num_bead_types,
  num_beads = num_beads
)
group2_beads <- create_attribute_matrix(
  c(
    rep((0.05 / 45), 45),
    rep((0.9 / 45), 45),  # 90% more likely to pick up certain beads than others
    rep((0.05 / 45), 45)
  ),
  num_sites = sample_size,
  num_bead_types = num_bead_types,
  num_beads = num_beads
)
group3_beads <- create_attribute_matrix(
  c(
    rep((0.05 / 45), 45),
    rep((0.05 / 45), 45),
    rep((0.9 / 45), 45) # 90% more likely to pick up certain beads than others
  ),
  num_sites = sample_size,
  num_bead_types = num_bead_types,
  num_beads = num_beads
)

atrribute_matrix <- rbind(group1_beads, group2_beads, group3_beads)
write.table(atrribute_matrix, "model2_binary_attributes.csv", row.names = FALSE, col.names = FALSE, sep = ",")
coordinate_matrix <- rbind(group1_coords, group2_coords, group3_coords)
write.table(coordinate_matrix, "model2_coordinates.csv", row.names = FALSE, col.names = FALSE, sep = ",")

# calculate disimilarity and geographical distances

jacquard_attributes <- as.matrix(vegan::vegdist(atrribute_matrix, method='jaccard'), binary = TRUE)
coordinate_distances <- great_arc_distance(coordinate_matrix)

# calculate and export mantel test results

sink("Model 2 Mantel Results.txt")
vegan::mantel(jacquard_attributes, coordinate_distances, permutations = 9999)
sink()

#

# model 3 - complete spatial and cultural seperation -------------

# set parameters

num_bead_types = 135
num_beads = 4
num_sites = 112
num_groups = 3
sample_size = num_sites / num_groups

# simulate spatial data

group1_coords <- simulate_regional_spatial_data(projected_europe, "Spain", crs, sample_size)
group2_coords <- simulate_regional_spatial_data(projected_europe, "Italy", crs, sample_size)
group3_coords <- simulate_regional_spatial_data(projected_europe, "Poland", crs, sample_size)

#visualise spatial data

plot(projected_europe, main = " ", asp = 1, type = "l")
points(group1_coords[,1], group1_coords[,2], pch = 19, cex = 1.5, col = "black")
points(group2_coords[,1], group2_coords[,2], pch = 19, cex = 1.5, col = "red")
points(group3_coords[,1], group3_coords[,2], pch = 19, cex = 1.5, col = "#4d72b2")

# simulate matrices for each group

group1_beads <- create_attribute_matrix(
  c(
    rep((1 / 45), 45),
    rep((0 / 45), 45), # probability of picking up a bead from another group is 0
    rep((0 / 45), 45) # probability of picking up a bead from another group is 0
  ),
  num_sites = sample_size,
  num_bead_types = num_bead_types,
  num_beads = num_beads
)
group2_beads <- create_attribute_matrix(
  c(
    rep((0 / 45), 45),
    rep((1 / 45), 45), # probability of picking up a bead from another group is 0
    rep((0 / 45), 45) # probability of picking up a bead from another group is 0
  ),
  num_sites = sample_size,
  num_bead_types = num_bead_types,
  num_beads = num_beads
)
group3_beads <- create_attribute_matrix(
  c(
    rep((0 / 45), 45),
    rep((0 / 45), 45), # probability of picking up a bead from another group is 0
    rep((1 / 45), 45) # probability of picking up a bead from another group is 0
  ),
  num_sites = sample_size,
  num_bead_types = num_bead_types,
  num_beads = num_beads
)

atrribute_matrix <- rbind(group1_beads, group2_beads, group3_beads)
write.table(atrribute_matrix, "model3_binary_attributes.csv", row.names = FALSE, col.names = FALSE, sep = ",")
coordinate_matrix <- rbind(group1_coords, group2_coords, group3_coords)
write.table(coordinate_matrix, "model3_coordinates.csv", row.names = FALSE, col.names = FALSE, sep = ",")

# calculate disimilarity and geographical distances

jacquard_attributes <- as.matrix(vegan::vegdist(atrribute_matrix, method='jaccard'), binary = TRUE)
coordinate_distances <- great_arc_distance(coordinate_matrix)

# calculate and export mantel test results

sink("Model 3 Mantel Results.txt")
vegan::mantel(jacquard_attributes, coordinate_distances, permutations = 9999)
sink()

#

# model 4 - Random spatial and non-random cultural attributes -------------

# set parameters

num_bead_types = 135
num_beads = 4
num_sites = 114
num_groups = 3
sample_size = num_sites / num_groups

# simulate spatial data

group1_coords <- simulate_random_spatial_data(projected_europe, sample_size)
group2_coords <- simulate_random_spatial_data(projected_europe, sample_size)
group3_coords <- simulate_random_spatial_data(projected_europe, sample_size)

#visualise spatial data

plot(projected_europe, main = " ", asp = 1, type = "l")
points(group1_coords[,1], group1_coords[,2], pch = 19, cex = 1.5, col = "black")
points(group2_coords[,1], group2_coords[,2], pch = 19, cex = 1.5, col = "red")
points(group3_coords[,1], group3_coords[,2], pch = 19, cex = 1.5, col = "#4d72b2")

# simulate matrices for each group

group1_beads <- create_attribute_matrix(
  c(
    rep((1 / 45), 45),
    rep((0 / 45), 45), # probability of picking up a bead from another group is 0
    rep((0 / 45), 45) # probability of picking up a bead from another group is 0
  ),
  num_sites = sample_size,
  num_bead_types = num_bead_types,
  num_beads = num_beads
)
group2_beads <- create_attribute_matrix(
  c(
    rep((0 / 45), 45),
    rep((1 / 45), 45), # probability of picking up a bead from another group is 0
    rep((0 / 45), 45) # probability of picking up a bead from another group is 0
  ),
  num_sites = sample_size,
  num_bead_types = num_bead_types,
  num_beads = num_beads
)
group3_beads <- create_attribute_matrix(
  c(
    rep((0 / 45), 45),
    rep((0 / 45), 45), # probability of picking up a bead from another group is 0
    rep((1 / 45), 45) # probability of picking up a bead from another group is 0
  ),
  num_sites = sample_size,
  num_bead_types = num_bead_types,
  num_beads = num_beads
)

atrribute_matrix <- rbind(group1_beads, group2_beads, group3_beads)
write.table(atrribute_matrix, "model4_binary_attributes.csv", row.names = FALSE, col.names = FALSE, sep = ",")
coordinate_matrix <- rbind(group1_coords, group2_coords, group3_coords)
write.table(coordinate_matrix, "model4_coordinates.csv", row.names = FALSE, col.names = FALSE, sep = ",")

# calculate disimilarity and geographical distances

jacquard_attributes <- as.matrix(vegan::vegdist(atrribute_matrix, method='jaccard'), binary = TRUE)
coordinate_distances <- great_arc_distance(coordinate_matrix)

# calculate and export mantel test results

sink("Model 4 Mantel Results.txt")
vegan::mantel(jacquard_attributes, coordinate_distances, permutations = 9999)
sink()

#

# model 5 - Complete Spatial Separation and random cultural attributes -------------

# set parameters

num_bead_types = 135
num_beads = 4
num_sites = 112
num_groups = 3
sample_size = num_sites / num_groups

# simulate spatial data

group1_coords <- simulate_regional_spatial_data(projected_europe, "Spain", crs, sample_size)
group2_coords <- simulate_regional_spatial_data(projected_europe, "Italy", crs, sample_size)
group3_coords <- simulate_regional_spatial_data(projected_europe, "Poland", crs, sample_size)

#visualise spatial data

plot(projected_europe, main = " ", asp = 1, type = "l")
points(group1_coords[,1], group1_coords[,2], pch = 19, cex = 1.5, col = "black")
points(group2_coords[,1], group2_coords[,2], pch = 19, cex = 1.5, col = "red")
points(group3_coords[,1], group3_coords[,2], pch = 19, cex = 1.5, col = "#4d72b2")

# simulate matrices for each group

group1_beads <- create_attribute_matrix(
  rep((1 / 135), 135), # the probability of drawing any type of bead is the same
  num_sites = sample_size,
  num_bead_types = num_bead_types,
  num_beads = num_beads
)
group2_beads <- create_attribute_matrix(
  rep((1 / 135), 135), # the probability of drawing any type of bead is the same
  num_sites = sample_size,
  num_bead_types = num_bead_types,
  num_beads = num_beads
)
group3_beads <- create_attribute_matrix(
  rep((1 / 135), 135), # the probability of drawing any type of bead is the same
  num_sites = sample_size,
  num_bead_types = num_bead_types,
  num_beads = num_beads
)

atrribute_matrix <- rbind(group1_beads, group2_beads, group3_beads)
write.table(atrribute_matrix, "model5_binary_attributes.csv", row.names = FALSE, col.names = FALSE, sep = ",")
coordinate_matrix <- rbind(group1_coords, group2_coords, group3_coords)
write.table(coordinate_matrix, "model5coordinates.csv", row.names = FALSE, col.names = FALSE, sep = ",")

# calculate disimilarity and geographical distances

jacquard_attributes <- as.matrix(vegan::vegdist(atrribute_matrix, method='jaccard'), binary = TRUE)
coordinate_distances <- great_arc_distance(coordinate_matrix)

# calculate and export mantel test results

sink("Model 5 Mantel Results.txt")
vegan::mantel(jacquard_attributes, coordinate_distances, permutations = 9999)
sink()

#
