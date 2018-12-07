library(lme4)
library(tidyverse)

# Import data
load('ldtsct.RData')

# Keep incorrect responses in data

# RT Exclusions
LDTSCT_DATA <- LDTSCT_DATA %>%
  filter(RT <= mean(LDTSCT_DATA$RT, na.rm=T)+2*sd(LDTSCT_DATA$RT, na.rm=T)) %>%
  filter(RT >= 250)

# Mean-centre and deviation-code task
LDTSCT_DATA$taskDev <- scale(ifelse(LDTSCT_DATA$task=="sct", 0, 1), scale = FALSE)

# Mean-centre and deviation-code stimuli (category-irrelevant words used as model's null)
LDTSCT_DATA$cwDev <- scale(ifelse(LDTSCT_DATA$stimType=="categ_word", 0, 1), scale = FALSE)  # category-relevant words
LDTSCT_DATA$pwDev <- scale(ifelse(LDTSCT_DATA$stimType=="pw_match", 0, 1), scale = FALSE)  # pseudowords
LDTSCT_DATA$nwDev <- scale(ifelse(LDTSCT_DATA$stimType=="nw_random", 0, 1), scale = FALSE)  # nonwords

# Mean-centre contrast-code stimuli for simple main effects (SMEs)
LDTSCT_DATA$cwCnt <- as.double(scale(ifelse(LDTSCT_DATA$stimType=="categ_word", 0, ifelse(LDTSCT_DATA$stimType=="nw_random", -0.5, 0.5)), scale = FALSE))  # category-relevant words
LDTSCT_DATA$pwCnt <- as.double(scale(ifelse(LDTSCT_DATA$stimType=="pw_match", 0, ifelse(LDTSCT_DATA$stimType=="categ_word", -0.5, 0.5)), scale = FALSE))  # pseudowords
LDTSCT_DATA$nwCnt <- as.double(scale(ifelse(LDTSCT_DATA$stimType=="nw_random", 0, ifelse(LDTSCT_DATA$stimType=="pw_match", -0.5, 0.5)), scale = FALSE))  # nonwords

# Mean-centre and deviation-code features of design
LDTSCT_DATA$respgroup <- scale(ifelse(LDTSCT_DATA$respgroup==1, 0, 1), scale = FALSE)  # 0 = right-as-yes, 1 = left-as-yes
LDTSCT_DATA$task_order <- as.double(scale(ifelse(LDTSCT_DATA$subgroup=="A"|LDTSCT_DATA$subgroup=="C", 0, 1), scale = FALSE))
LDTSCT_DATA$stim_list <- as.double(scale(ifelse(LDTSCT_DATA$subgroup=="A"|LDTSCT_DATA$subgroup=="B", 0, 1), scale = FALSE))

# Full model formula
form <- corr_ans ~ taskDev * (cwDev + nwDev + pwDev) + task_order + stim_list + respgroup +
  (1 + taskDev * (cwDev + nwDev + pwDev) || subid) +
  (1 + taskDev || item)

# Write output to file
sink(file = sprintf("results_%s.txt", as.character(form[2])), type = c("output"))

# Fit full model
m.full <- glmer(form, data = LDTSCT_DATA, family = binomial(logit), control = glmerControl(optimizer = c("bobyqa")))

# Print summary of model
summary(m.full)

# Main Effect of task
m.task <- update(m.full, .~. -taskDev)
anova(m.full, m.task, test="Chi")

# Main Effect of stimulus type
m.stim <- update(m.full, .~. -cwDev -nwDev -pwDev)
anova(m.full, m.stim, test="Chi")

# Effect of category-relevant words
m.stim.cw <- update(m.full, .~. -cwDev)
anova(m.full, m.stim.cw, test="Chi")

# Effect of pseudowords
m.stim.pw <- update(m.full, .~. -pwDev)
anova(m.full, m.stim.pw, test="Chi")

# Effect of nonwords
m.stim.nw <- update(m.full, .~. -nwDev)
anova(m.full, m.stim.nw, test="Chi")

# Task-stimulus interaction
m.int <- update(m.full, .~. -taskDev:cwDev -taskDev:pwDev -taskDev:nwDev)
anova(m.full, m.int, test="Chi")

################################################################
# Effect of design

# Effect of task order
m.to <- update(m.full, .~. -task_order)
anova(m.full, m.to, test="Chi")

# Effect of stimulus list
m.sl <- update(m.full, .~. -stim_list)
anova(m.full, m.sl, test="Chi")

# Effect of response group
m.rg <- update(m.full, .~. -respgroup)
anova(m.full, m.rg, test="Chi")

################################################################
# SME - Effect of task at each level of stimulus

# Dummy-code stimuli for simple main effects (SMEs)
LDTSCT_DATA$cwDum <- ifelse(LDTSCT_DATA$stimType=="categ_word", 1, 0)  # category-relevant words
LDTSCT_DATA$pwDum <- ifelse(LDTSCT_DATA$stimType=="pw_match", 1, 0)  # pseudowords
LDTSCT_DATA$nwDum <- ifelse(LDTSCT_DATA$stimType=="nw_random", 1, 0)  # nonwords

# Get formulae for dummy-coded stimuli

eval(parse(text=sprintf('form.cw <- %s ~ taskDev * cwDum + task_order + stim_list + respgroup +
                        (1 + taskDev * cwDum || subid) +
                        (1 + taskDev || item)
                        ', as.character(form[2]))))

eval(parse(text=sprintf('form.pw <- %s ~ taskDev * pwDum + task_order + stim_list + respgroup +
                        (1 + taskDev * pwDum || subid) +
                        (1 + taskDev || item)
                        ', as.character(form[2]))))

eval(parse(text=sprintf('form.nw <- %s ~ taskDev * nwDum + task_order + stim_list + respgroup +
                        (1 + taskDev * nwDum || subid) +
                        (1 + taskDev || item)
                        ', as.character(form[2]))))

# Examine effect of task for category-relevant words
m.cw <- glmer(form.cw, data = LDTSCT_DATA, family = binomial(logit), control = glmerControl(optimizer = c("bobyqa")))
m.cw.int <- update(m.cw, .~. -taskDev)
anova(m.cw, m.cw.int, test="Chi")

# Examine effect of task for pseudowords
m.pw <- glmer(form.pw, data = LDTSCT_DATA, family = binomial(logit), control = glmerControl(optimizer = c("bobyqa")))
m.pw.int <- update(m.pw, .~. -taskDev)
anova(m.pw, m.pw.int, test="Chi")

# Examine effect of task for nonwords
m.nw <- glmer(form.nw, data = LDTSCT_DATA, family = binomial(logit), control = glmerControl(optimizer = c("bobyqa")))
m.nw.int <- update(m.nw, .~. -taskDev)
anova(m.nw, m.nw.int, test="Chi")


# Stop printing output to file
sink()
