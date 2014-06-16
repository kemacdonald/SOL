
iChart <- iChart[!is.na(iChart$RT),]

for(i_row in 1:nrow(iChart)){
	
	for(i_col in 147:207){
		
		if(as.numeric(names(iChart)[i_col]) < as.numeric(iChart$RT[i_row])){
			
			iChart[i_row, i_col] <- 0
			
			}else{
			
			iChart[i_row, i_col] <- 1
			
			 }
			iChart[,i_col] <- as.numeric(iChart[,i_col])
		}	
		print(i_row)
	}
	#207