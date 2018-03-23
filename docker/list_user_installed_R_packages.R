# J. Taroni 2018
# Get list of installed R packages on Docker image
# Adapted from:
# https://www.r-bloggers.com/list-of-user-installed-r-packages-and-their-versions/
# https://stackoverflow.com/questions/38481980/get-the-list-of-installed-packages-by-user-in-r

docker.tag <- "recount"
dir.create(file.path("docker", docker.tag), recursive = TRUE,
           showWarnings = FALSE)

packages.df <- as.data.frame(installed.packages()[, c(1, 3:4)])
packages.df <- packages.df[is.na(packages.df$Priority), 1:2, drop=FALSE]

sink(file.path("docker", docker.tag, "user_installed_R_packages.txt"))
print(packages.df, row.names=FALSE)
sink()
