genes = read.table("batch2.output.genes.txt", header=T)


# count how many genes are included in the analysis
dim(genes) 

# take a look at the first few rows in the “genes” matrix:
head(genes)

meanTPM = data.frame(ID=genes[,1], Means=rowMeans(genes[,-1]))

genes = genes[meanTPM[,2]>0, ]  ## only keep the genes that have non-zero expression in at least one sample

## exploratory figures: 
# histogram of gene expression values (TPM)

hist(log10(genes[,2]  + 1e-4), xlab="log10(gene TPM)")

# mark the median value
abline(v=log10(median((genes[,2]))))


# utility functions to calculate Jensen-Shannon divergence for doing multi-dimensional scale (MDS) analysis: 

shannon.entropy <- function(p)
{
        if (min(p) < 0 || sum(p) <= 0)
                return(NA)
        p.norm <- p[p>0]/sum(p)
        -sum(log2(p.norm)*p.norm)
}

jensen_shannon <- function(p, q){
        ## JSD = H(0.5 *(p +q)) - 0.5H(p) - 0.5H(q)
        # H(X) = \sum x_i * log2(x_i)

        p = p / sum(p)
        q = q / sum(q)
        Hj = shannon.entropy(0.5 *(p+q)) 
        Hp = shannon.entropy(p) 
        Hq = shannon.entropy(q)
        
        jsd = Hj - 0.5*(Hp+Hq)
#       cat(Hj, Hp, Hq, jsd, "\n")
        return(jsd)
}

## take genes with reasonble expression
# y  = genes[genes$x5 > 0.5, ]


y = genes

## MDS plot based on Jensen-Shannon Metric
nsample = ncol(y) - 1
jsd = matrix(nrow=nsample, ncol=nsample)
colnames(jsd) = colnames(y)[2: (nsample + 1)]
rownames(jsd) = colnames(y)[2: (nsample + 1)]

for (i in 1: (nsample -1) ) {
	jsd[i,i] = 0
	for (j in (i+1):nsample) {
		jsd[j,j] = 0
	##	print(paste(i, j, sep=" "))
		jsd[i,j] = jensen_shannon(y[,i+1], y[,j+1])
		jsd[j,i] = jsd[i,j]
	}
}

## square root of JSD is JSM, which is a metric that can be used as distance. 
# do MDS with 2 dimensions
mdsfit <- cmdscale(jsd**0.5, eig = T, k = 2 )

plot(mdsfit$points,  xlab="MDS D1", ylab="MDS D2", pch=21, col='blue', lwd=2)

grid()

text(mdsfit$points, label=rownames(mdsfit$points), adj = -0.2, col="gray")

### the MDS plot is attached in the email. Samples 5 and 6 are controls. Dimension 1 separates Samples 3,4,5,6 from other rest, which may make sense.  I’m not sure how to interpret the result of sample 1 and 2. 

# sanity check of key genes:

genes[genes$geneName == "Ctnnb1", ]
genes[genes$geneName == "Ep300", ]


## simple hierarchical clustering: 
d = dist(t(y[,2: (nsample + 1)]), )
h = hclust(d)
plot(h)


## revisit the main cluster: 
points = mdsfit$points

mainc = points[points[,1] < 0.1 & points[,2] < 0.1, ]
maincSamples = rownames(mainc)

newjsd = jsd[maincSamples, ]
newjsd = newjsd[,maincSamples]

newmdsfit <- cmdscale(newjsd**0.5, eig = T, k = 2 )

plot(newmdsfit$points, xlab="MDS D1", ylab="MDS D2", pch=21, col='red', lwd=2, main="main cluster")
text(newmdsfit$points, label=rownames(newmdsfit$points), adj = -0.2, col="gray")

