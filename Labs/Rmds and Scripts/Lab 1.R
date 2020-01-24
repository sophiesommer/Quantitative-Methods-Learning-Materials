##################################################
# Intro to R and T Tests, Simple Linear Regression
##################################################

# A few basics:

# To "comment out" text, use the # symbol or highlight chunks of text and use shift + command + C
# To run a line of code, use command + enter
# To get more information, use a ? before a command and run it (e.g. ?setwd) or go to the Help tab in the lower right window

# Text wrapping: Tools > Global Options > Code (Editing tab) > Soft-wrap R source files
# Change colors: Tools > Global Options > Appearance

#########################################
# 1. Get working directory and set working directory
# Your working directory is the folder you are working in and where new files will be created and saved
getwd()
setwd("/Users/sophiesommer/Desktop/Grad School/2003 TA Materials")
# *** A shortcut to finding the file path is to go to the folder in Finder, right-click and choose "Get Info", and then copy and paste the file path from "Where:" under "General" ***

#########################################
# 2. Create, open, or save an R script
## CREATE
# - File > New File > R Script
# - New file icon at top left
# - shift + command + N

## OPEN
# - File > Open File...
# - command + O

## SAVE
# - File > Save
# - command + S


#########################################
# 3. Install a package and load a library
## INSTALL A PACKAGE
# Make sure that the package name is in quotes. You will also need internet access to download packages
install.packages("openintro")

# Load a library
library(openintro)

# OR: use require() which will first check if the library is loaded before re-loading it
require(openintro)

#########################################
# 4. Load a data set and inspect
## LOADING DATA
# read.csv() loads csv (comma-separated-values) files (e.g. Excel files)
# Here, I've named the dataset marathon.csv and saved it to my Data folder
marathon <- read.csv("./Data/marathon.csv")

# The next option is great if you're feeling lazy or desperate, but not so great for running R markdown files or anything you want someone else to be able to run...
marathon <- read.csv(file.choose())

# There are also some built in datasets in R (some are also available via packages). 
# The marathon dataset above can also be accessed via the openintro package. So, once you've run the code above to load the library, you can load the dataset simply by using the data() function:
data(marathon)

# To delete an object use rm(). BE VERY CAREFUL WITH DELETING THINGS!
rm(marathon)


## INSPECTING DATA
# The follwing command is the same as clicking the dataframe icon in the Environment window (upper right)
View(marathon)

# Check what type of object income is (certain commands only work for certain classes of objects). This is generally not necessary to check every time you open a dataset; it's just for future reference**
class(marathon)

# Inspect the first and last rows of the dataset
head(marathon) # by default, shows first 6 rows
tail(marathon, n = 2) # tail shows last rows, n=2 means the last 2 rows

#########################################
# 5. Check some descriptives of the variables
# To reference a specific column use dataframe_name$column

# Just the first 6 values of the Gender variable
head(marathon$Gender)
# Counts of each value of Gender
table(marathon$Gender)
# Lowest and highest values of the Time variable
range(marathon$Time)

# Quick summary stats of variables
summary(marathon$Time)
summary(marathon$Time, digits = 2)

# Calculate individual stats
mean(marathon$Time) #mean
var(marathon$Time) #variance
sd(marathon$Time) #standard deviation
sqrt(var(marathon$Time)) #note that sd=sqrt(var)
median(marathon$Time) #median
cor(marathon$Year, marathon$Time) #correlation
cov(marathon$Year, marathon$Time) #covariance
# Question: How could you interpret the negative correlation/covariance between year and marathon time?

#########################################
# 6. A trick to avoid typing marathon$ before variable names every time

# Two options:  
# 1. We can use the attach() function to "attach" the dataset to the file path that R searches when you reference variables.
# Example:
attach(marathon)
mean(Time)
# To check which dataframes are attached:
search()
# To detach a dataset
detach(marathon)
mean(Time)
# Attach is a useful tool, but the issue here is that it can get confusing if you're working with more than one dataset (especially if those datasets have similar variable names), or if you forget to detach the dataset.

# 2. Save them as separate variables
Time <- marathon$Time
Gender <- marathon$Gender

#########################################
# 7. Setting up for a T-test (or Simple Regression)
# For this part, we will run a simple linear regression of time as a function of gender.
# In other words, we want to know whether marathon times differ by gender

## First, let's check the class of the variables that we will use. 
class(marathon$Gender)
class(marathon$Time)
# Note that gender is a factor variable with two levels, f and m. You can also use levels to see the order:
levels(marathon$Gender)
#f comes first, so essentially, f=0 and m=1 when we run the linear regression

# Note: If you want to make this more explicit, we will change all of the f's to 0s and m's to 1's 
# (this is not necessary in order to run the model, but will be helpful for plotting):
marathon$Gender <- ifelse(marathon$Gender=="f", 0, 1)

# Re-inspect the gender variable to confirm it is now numeric with 0s and 1s
class(marathon$Gender) #inspect class
table(marathon$Gender) #look at a table of values: there are 29 0s (females) and 30 1s (males)

# Subset the data into two vectors: marathon times for females and marathon times for males
female_times = marathon$Time[marathon$Gender==0]
male_times <- marathon$Time[marathon$Gender==1]

## Calculate N, mean, variance, and SE for each group

# N
female_N <- length(female_times) 
female_N
male_N <- length(male_times)
male_N

# Mean
female_mean <- mean(female_times) 
female_mean
male_mean <- mean(male_times)
male_mean

# Variance
female_variance <- var(female_times) 
female_variance
male_variance <- var(male_times)
male_variance

# SE
female_SE <- sqrt(female_variance/female_N)
female_SE
male_SE <- sqrt(male_variance/male_N)
male_SE

## Levene test
install.packages("car")
require(car)
# Note: need to use factor() to coerce Gender back to a factor variable for purposes of this test
# Levene test tests H0: equal variance, HA: non-equal variance
leveneTest(Time~factor(Gender), data=marathon)
# p=.1589 so at a .05 significance level, we cannot reject the null that the variances are equal
# so we can use equal variances for t test

# T test assuming equal variances
# This tests H0: mean female time= mean male time vs HA: mean female time not equal to mean male time
t.test(Time~Gender, data=marathon, var.equal=TRUE)
# alternatively (does the same thing):
t.test(female_times, male_times, var.equal=TRUE)

# If you wanted to run the T test not assuming equal variances:
t.test(Time~Gender, data=marathon, var.equal=FALSE)

#########################################
# 8. Running a simple linear regression  

# There are LOTS of arguments you can specify, but we just need the basics, so most of the default setting should be fine.
lin.model <- lm(Time~Gender, data=marathon)
# First we specifiy the model. Outcome/response variable goes on the left, explanatory variables go on the right.
# We don't need to reference the dataset and use $'s because we've specified the dataset with the argument data=marathon.
# We can verify that this object is a list using is.list()
is.list(lin.model)
# We can see all the first level items in the list using ls()
ls(lin.model)
lin.model$coefficients
# Using just "coef" as done in the lecture notes is simply a shortcut for typing out the whole word

summary(lin.model)
# You'll be using summary() on models more as you progress and need to evaluate how well your model fits the data or when you might decide to go with one model over another

plot(lin.model)
# You'll be able to use the plot() function on your model to check assumptions about that model, such as whether it meets assumptions of normality or homoscedasticity. Again, more on this later!

#########################################
# 9. Plotting your data and your model
# Usually we want to plot our data early on when doing descriptives so we have a visual of what's going on in the data. We'll do a quick scatter plot here though.
plot(marathon$Gender, marathon$Time)

# We can add better labels to the plot and a title
plot(marathon$Gender, marathon$Time, main = "Simple Linear Regression Example", xlab = "Gender", ylab = "Time")
# We can also add our newly found regression line
abline(lin.model, col = "green")

# We can also add in the points (0, mean female time) and (1, mean male time)
# The regression line will run through these points
points(0, female_mean, col = 6, pch = 19)
points(1, male_mean, col = 6, pch = 19)

