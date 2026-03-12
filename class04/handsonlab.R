#### vectors ####

## types of vectors

# numeric
x <- 1:6

# character
y <- c("alice", "chandra", "amy", "barry")
paste(y, "loves R") # adds 2nd argument to each element

# logical
z <- c(TRUE, FALSE, TRUE, TRUE, FALSE)
z + 100 # works because TRUE = 1, FALSE = 
x > 3 # returns F F F T T T
y == "alice" # returns T F F F

#### data frames ####

# creating a dataframe from scratch
df <- data.frame(
  nums = 1:5,
  chars = letters[1:5],
  logical = c(T, T, F, T, F))

# selecting some parts of dataframe
df[3, ] # prints out entire 3rd row
df[, 2] # prints out entire 2nd column

df$chars # alternative to above
df$chars[3] # selects element in 3rd row of chars column
