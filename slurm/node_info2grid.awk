BEGIN {
    need_header=1
}
# parse data
($1!=""){
    split($0,nvp," ");
    for (x in nvp) { 
	#print $x; 
	split($x,nv,"="); 
	#print nv[1] ">>" nv[2]; 
	vals[nv[1]]=nv[2]
	if( need_header > 0) {
	    names[need_header] = nv[1]
	    need_header++
	}
n    } 
}
# EMIT column headers
($1=="" && need_header>0){
#    printf "#" 
    #for(x in vals){printf x"\t"};
    for(i =1; i<=length(vals);i++) { printf names[i] "\t"};
    printf "\n" ;
    need_header=0;
}
# EMIT values for a node in a row
($1=="" && vals["NodeName"] != ""){
    #for(x in vals){printf vals[x]"\t" }; 
    for(i =1; i<=length(vals);i++) { printf vals[names[i]] "\t"};
    printf "\n"
    delete vals
}
