###

### assembleComm

### X=environmental parameter data (selected from offering)
### fun=vector of selected species functions


assembleComm <- function(x, fun){

saa<-matrix(,nrow=3, ncol=2) #build matrix

for(i in 1:length(fun)) {
  	for(row in 1:nrow(df)){
   saa[row,i]<-do.call(fun[[i]], list(df[row,1],df[row,2]))
      }
	}



}
