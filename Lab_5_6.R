# Czyszczenie od anomalii


library(rgl)
library(scatterplot3d)
library(dplyr)
library(ggplot2)



set.seed(123)  # dla powtarzalności
df1 = read.csv("khomenko_output_bez_x.txt", sep=";", dec=".")

str(df1)
print(min(df1$Z1))
print(min(df1$Z2))
print(min(df1$Z3))

print(max(df1$Z1))
print(max(df1$Z2))
print(max(df1$Z3))




# ponowne usuwanie anomalii
# Odcinamy po 2.5% z każdej strony
xmin <- quantile(df1$X_C, 0.025)
xmax <- quantile(df1$X_C, 0.975)
ymin <- quantile(df1$Y_C, 0.025)
ymax <- quantile(df1$Y_C, 0.975)

# Zachowujemy punkty wewnętrzne (czyli środek)
df_trimmed <- subset(df1,
                     X_C >= xmin & X_C <= xmax &
                       Y_C >= ymin & Y_C <= ymax)
print(nrow(df_trimmed))
nrow(df)






grupowanie_normals_N<-kmeans(as.matrix(dplyr::select(df_trimmed, c("X_N","Y_N", "Z_N"))), centers
                           = 4, nstart=40, iter.max = 100000, algorithm="Lloyd" )

str(grupowanie_normals_N)

norm_vec <- function(x) sqrt(sum(x^2))


transform_func_N <- function(centers_matrix) {
  result <- data.frame(
    Dip_dir = numeric(),
    Dip_ang = numeric(),
    clustering = character()
  )
  for (i in 1:nrow(centers_matrix)) {
    row_curr = centers_matrix[i, ]
    row_curr <- row_curr/norm_vec(row_curr)
    x <- row_curr[1]
    y <- row_curr[2]
    z <- row_curr[3]
    
    dip_dir <- atan2(y, x) * 180/pi
    if (dip_dir < 0) dip_dir <- (dip_dir + 360)%%360
    else dip_dir <- dip_dir%%360
    dip_ang = acos(z) * 180/pi
    result <- rbind(result, data.frame(
      Dip_dir = dip_dir,
      Dip_ang = dip_ang,
      clustering = i
    ))
  }
  return(result)
}



centers_N = grupowanie_normals_N$centers
centers_N_calc = transform_func_N(centers_N)
centers_N_calc
write.table(centers_N_calc[1,],row.names=FALSE,col.names=TRUE,sep=",", quote=FALSE,
            file="01_N.txt")
write.table(centers_N_calc[2,],row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE,
            file="02_N.txt")
write.table(centers_N_calc[3,],row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE,
            file="03_N.txt")
write.table(centers_N_calc[4,],row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE,
            file="04_N.txt")


df_trimmed$clustering_N <- grupowanie_normals_N$cluster
df_trimmed

str(df_trimmed)
unique(df_trimmed)
table(df_trimmed$clustering, useNA = "ifany")



# dla 20000
sample_df = sample_n(df_trimmed, 5000)


df1 <- sample_df[sample_df$clustering_N == 1,]
df1
df2 <- sample_df[sample_df$clustering_N == 2,]
df3 <- sample_df[sample_df$clustering_N == 3,]
df4 <- sample_df[sample_df$clustering_N == 4,]

nrow(df1)
nrow(df2)
nrow(df3)
nrow(df4)

write.table(df1, file = "cluster_N_1.txt", sep = ";", dec = ".", row.names = FALSE, quote = FALSE)
write.table(df2, file = "cluster_N_2.txt", sep = ";", dec = ".", row.names = FALSE, quote = FALSE)
write.table(df3, file = "cluster_N_3.txt", sep = ";", dec = ".", row.names = FALSE, quote = FALSE)
write.table(df4, file = "cluster_N_4.txt", sep = ";", dec = ".", row.names = FALSE, quote = FALSE)


ggplot(df_trimmed, aes(x = Y_C, y = X_C, color = factor(clustering_N))) +
  geom_point(size = 0.3) +
  labs(color = "Klastry orientacji powierzchni")+
  coord_equal() +
  scale_color_manual(
    values = c(
      "1" = "#a1d99b", # zielone
      "2" = "#fdae6b", # pomaranczowy
      "3" = "#2b8cbe", # niebieski
      "4" = "#756bb1"  # fioletowy
    )
  ) + 
  labs(title = "Grupowanie wektorów normalnych (k-means)",
       x = "X_C", y = "Y_C") +
  theme_minimal()+ coord_equal()

table(df_trimmed$clustering, useNA = "ifany")


#df3["Dip_ang"]

stereonet_df1 = df1[, c("Dip_dir", "Dip_ang")]
stereonet_df1
nrow(stereonet_df1)
write.table(stereonet_df1, file = "stereonet_cluster_N_1_subset_5000.txt", sep = ",", dec = ".", row.names = FALSE, quote = FALSE)

stereonet_df2 = df2[, c("Dip_dir", "Dip_ang")]
write.table(stereonet_df2, file = "stereonet_cluster_N_2_subset_5000.txt", sep = ",", dec = ".", row.names = FALSE, quote = FALSE)

stereonet_df3 = df3[, c("Dip_dir", "Dip_ang")]
write.table(stereonet_df3, file = "stereonet_cluster_N_3_subset_5000.txt", sep = ",", dec = ".", row.names = FALSE, quote = FALSE)

stereonet_df4 = df4[, c("Dip_dir", "Dip_ang")]
write.table(stereonet_df4, file = "stereonet_cluster_N_4_subset_5000.txt", sep = ",", dec = ".", row.names = FALSE, quote = FALSE)


nrow(stereonet_df1)
nrow(stereonet_df2)
nrow(stereonet_df3)
nrow(stereonet_df4)
#nrow(stereonet_df3)


stereonet_df2
hist(stereonet_df2$Dip_dir)





# Dla X_D, Y_D, Z_D
grupowanie_normals_D<-kmeans(as.matrix(dplyr::select(df_trimmed, c("X_D","Y_D", "Z_D"))), centers
                             = 4, nstart=40, iter.max = 100000, algorithm="Lloyd" )

str(grupowanie_normals_D)

norm_vec <- function(x) sqrt(sum(x^2))


transform_func_D <- function(centers_matrix) {
  result <- data.frame(
    Dip_dir = numeric(),
    Dip_ang = numeric(),
    clustering = character()
  )
  for (i in 1:nrow(centers_matrix)) {
    row_curr = centers_matrix[i, ]
    row_curr <- row_curr/norm_vec(row_curr)
    x <- row_curr[1]
    y <- row_curr[2]
    z <- row_curr[3]
    
    dip_dir <- atan2(y, x) * 180/pi
    if (dip_dir < 0) dip_dir <- (dip_dir + 360)%%360
    else dip_dir <- dip_dir%%360
    dip_ang = acos(z) * 180/pi-90
    result <- rbind(result, data.frame(
      Dip_dir = dip_dir,
      Dip_ang = dip_ang,
      clustering = i
    ))
  }
  return(result)
}


# wynik = transform_func_D(a)
# wynik

centers_D = grupowanie_normals_D$centers
centers_D_calc = transform_func_D(centers_D)
centers_D_calc
write.table(centers_D_calc[1,],row.names=FALSE,col.names=TRUE,sep=",", quote=FALSE,
            file="01_D.txt")
write.table(centers_D_calc[2,],row.names=FALSE,col.names=TRUE,sep=",", quote=FALSE,
            file="02_D.txt")
write.table(centers_D_calc[3,],row.names=FALSE,col.names=TRUE,sep=",", quote=FALSE,
            file="03_D.txt")
write.table(centers_D_calc[4,],row.names=FALSE,col.names=TRUE,sep=",", quote=FALSE,
            file="04_D.txt")


df_trimmed$clustering_D <- grupowanie_normals_D$cluster
df_trimmed

str(df_trimmed)
unique(df_trimmed)
table(df_trimmed$clustering, useNA = "ifany")



# dla 5000
sample_df = sample_n(df_trimmed, 5000)


df1 <- sample_df[sample_df$clustering_D == 1,]
df1
df2 <- sample_df[sample_df$clustering_D == 2,]
df3 <- sample_df[sample_df$clustering_D == 3,]
df4 <- sample_df[sample_df$clustering_D == 4,]

nrow(df1)
nrow(df2)
nrow(df3)
nrow(df4)

write.table(df1, file = "cluster_D_1.txt", sep = ";", dec = ".", row.names = FALSE, quote = FALSE)
write.table(df2, file = "cluster_D_2.txt", sep = ";", dec = ".", row.names = FALSE, quote = FALSE)
write.table(df3, file = "cluster_D_3.txt", sep = ";", dec = ".", row.names = FALSE, quote = FALSE)
write.table(df4, file = "cluster_D_4.txt", sep = ";", dec = ".", row.names = FALSE, quote = FALSE)


ggplot(df_trimmed, aes(x = Y_C, y = X_C, color = factor(clustering_D))) +
  geom_point(size = 0.3) +
  labs(color = "Klastry kierunku spadku")+
  coord_equal() +
  scale_color_manual(
    values = c(
      "1" = "#a1d99b", # zielony
      "2" = "#fdae6b", # pomaranczowy
      "3" = "#2b8cbe", # niebieski
      "4" = "#756bb1"  # fioletowy
    )
  ) + 
  labs(title = "Grupowanie wektorów kierunku spadku (k-means)",
       x = "X_C", y = "Y_C") +
  theme_minimal()+ coord_equal()

table(df_trimmed$clustering_D, useNA = "ifany")


#df3["Dip_ang"]

stereonet_df1 = df1[, c("Dip_dir", "Dip_ang")]
nrow(stereonet_df1)
write.table(stereonet_df1, file = "stereonet_cluster_D_1_subset_5000.txt", sep = ",", dec = ".", row.names = FALSE, quote = FALSE)

stereonet_df2 = df2[, c("Dip_dir", "Dip_ang")]
write.table(stereonet_df2, file = "stereonet_cluster_D_2_subset_5000.txt", sep = ",", dec = ".", row.names = FALSE, quote = FALSE)

stereonet_df3 = df3[, c("Dip_dir", "Dip_ang")]
write.table(stereonet_df3, file = "stereonet_cluster_D_3_subset_5000.txt", sep = ",", dec = ".", row.names = FALSE, quote = FALSE)

stereonet_df4 = df4[, c("Dip_dir", "Dip_ang")]
write.table(stereonet_df4, file = "stereonet_cluster_D_4_subset_5000.txt", sep = ",", dec = ".", row.names = FALSE, quote = FALSE)


nrow(stereonet_df1)
nrow(stereonet_df2)
#nrow(stereonet_df3)


stereonet_df2
hist(stereonet_df2$Dip_dir)


