# read in the files
genome = read.table('mafs.txt', skip=1, header=TRUE)
annotation <- scan("mafs.txt", what="character", nlines=1, sep='!')
mutants = read.table('muts.txt', header=TRUE)

# make the figure
pdf('mutation_simulation.pdf',8,32)
op <- par(mfrow=c(4,1))
plot.ts(genome[,-1:-2], plot.type="single", col=sample(rainbow(10)), xlab="Generation", ylab="Minor Allele Frequency")
mtext(annotation, side=3)
barplot(genome$Pop, new=T, names.arg=genome$Gen, xlab="Generation", ylab="Population Size")
boxplot(Mut~Gen, data=mutants, col="lightblue", xlab="Generation", ylab="Mutation Frequency")
boxplot(log10(Fit)~Gen, data=mutants, col="mistyrose", xlab="Generation", ylab="Fitness")
abline(a=0, b=0)
dev.off()

# zeroes <- subset(mutants, Mut==0)
# table(zeroes)


#library(reshape2)
#new <- melt(genome, id.vars=c("Gen", "Pop"), measure.vars=3:dim(genome)[2], variable.name="Pos", value.name="MAF")

#ggplot(new, aes(x=Gen, y=MAF, colour=Pos)) +
#  geom_line() +
#  scale_x_continuous(breaks=0:dim(genome)[1]) +
#  scale_fill_brewer(guide=FALSE) +
#  geom_bar()