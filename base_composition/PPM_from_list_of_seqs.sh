
cat ${fasta} |
  grep -v '>' |
  awk '{L=length($1);for(i=1;i<=L;i++) {B=substr($1,i,1);T[i][B]++;}} END{for(BI=0;BI<4;BI++) {B=(BI==0?"A":(BI==1?"C":(BI==2?"G":"T")));printf("%s",B); for(i in T) {tot=0.0;for(B2 in T[i]){tot+=T[i][B2];}printf("\t%0.2f",(T[i][B]/tot));} printf("\n");}}' 
