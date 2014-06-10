files <- list.files(pattern="\\.rds$", recursive=TRUE)

for (file in files)
{
    cat(paste("**", file, "**\n"))
    print(readRDS(file))
    cat("\n")
}
