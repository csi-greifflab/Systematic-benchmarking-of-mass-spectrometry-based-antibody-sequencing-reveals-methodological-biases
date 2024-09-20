# Read the list of filenames from list_of_files.txt
filelist_PGDM1400 <- readLines("./filelist_PGDM1400.txt")
filelist_Umab <- readLines("./filelist_Umab.txt")
filelist_h9C12WT <- readLines("./filelist_h9C12WT.txt")

# The original block of code as a string
original_code <- "Casanovo ->
        Path       : /storage/khangl/MS_benchmarking/Casanovo_output/EH3420.mztab
        Name       : EH3420
        CutoffScore: 0.8
        -FilterPPM  : 100
    <-"

# Create a function to replace "EH3420" with a given filename
generate_code_block <- function(filename) {
  new_code <- gsub("EH3420", filename, original_code)
  return(new_code)
}

# PGDM1400
# Initialize an empty string to accumulate all code blocks
all_code_blocks <- ""

# Loop through each filename and generate the new blocks of code
for (name in filelist_PGDM1400) {
  new_code_block <- generate_code_block(name)
  
  # Accumulate the new block of code into the single string
  all_code_blocks <- paste(all_code_blocks, new_code_block, "\n\n", sep = "")
}

# Define the output file name
output_filename <- "all_blocks_PGDM1400.txt"

# Write the accumulated code blocks to the single file
writeLines(all_code_blocks, output_filename)

cat("All code blocks written to", output_filename, "\n")

# Umab
# Initialize an empty string to accumulate all code blocks
all_code_blocks <- ""

# Loop through each filename and generate the new blocks of code
for (name in filelist_Umab) {
  new_code_block <- generate_code_block(name)
  
  # Accumulate the new block of code into the single string
  all_code_blocks <- paste(all_code_blocks, new_code_block, "\n\n", sep = "")
}

# Define the output file name
output_filename <- "all_blocks_Umab.txt"

# Write the accumulated code blocks to the single file
writeLines(all_code_blocks, output_filename)

cat("All code blocks written to", output_filename, "\n")

# h9C12WT
# Initialize an empty string to accumulate all code blocks
all_code_blocks <- ""

# Loop through each filename and generate the new blocks of code
for (name in filelist_h9C12WT) {
  new_code_block <- generate_code_block(name)
  
  # Accumulate the new block of code into the single string
  all_code_blocks <- paste(all_code_blocks, new_code_block, "\n\n", sep = "")
}

# Define the output file name
output_filename <- "all_blocks_h9C12WT.txt"

# Write the accumulated code blocks to the single file
writeLines(all_code_blocks, output_filename)

cat("All code blocks written to", output_filename, "\n")