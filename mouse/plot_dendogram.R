setwd('C:/Users/Rolf/Desktop/Dropbox/work connectomics/eclipse/mouse')
library(ape)
data=read.csv('parcellation_hierarchy.csv')
data=data[-which((data[,2])==0),]
n=dim(data)[1]
rootID=997

ancestor=matrix(0,n,1)#gives the direct ancestor for each node
ids=as.numeric(data[,1])
for(i in 1:n){
	strs=strsplit(as.character(data[i,3]),'/')[[1]];
if(length(strs)==2){
	id=-1;
}else
	id= which(ids==as.numeric(strs[length(strs)-1]))
if(length(id)==1){
	ancestor[i] = id
}
else
	if(rootID==as.numeric(strs[length(strs)-1])){
		ancestor[i]==-1;
	}
	else
		ancestor[i]=-2;
}




#try to make a phylo object
#a=ancestor
n=n-1
#root=which(ancestor[,1]==-1);#remove root
#a=a[-root]
#removed=root
#for(i in 1:(n-1)){
#	if(a[i]==root){
#		a[i]=n
#}
#	else if(a[i]>root)
#		a[i]=a[i]-1;
#}
remove=c()##purge all the singles
for(i in 1:n){
	if(sum(ancestor==ancestor[i])==1){
		if(ancestor[i]!=-1){#root
#			its offspring
			
##			add a dummy node, to preserve all info
#			data=rbind(data,c(-1,-1,-1,'Dummy','DUM',1,1,0,NA))
#			n=dim(data)[1]
#			ancestor[n]=ancestor[i];
			remove=c(remove,i)
		}
		
	}
#	remove layes (ie things with numbers)
		if(grepl('[0123456789]',data$Abbreviation[i])){
			remove=c(remove,i)
		}
}
#data$Volume..in.parcellation.=as.numeric(data$Volume..in.parcellation.)
#data$Volume..original.=as.numeric(data$Volume..original.)
#d=data;


#replace dorsal striatum by CP
remove[remove==603]=230

for(i in remove){
	parentid=ancestor[i]
	while(parentid%in%remove)
		parentid=ancestor[parentid]
	data$Volume..original.[parentid]=data$Volume..original.[parentid]+data$Volume..original.[i]
	
}
d=data[-remove,]
n=dim(d)[1]
ancestor=matrix(0,n,1)#gives the direct ancestor for each node
ids=as.numeric(d[,1])
for(i in 1:n){
	strs=strsplit(as.character(d[i,3]),'/')[[1]];
	if(length(strs)<=3){
		id=-1;
	}else{
	if(rootID==as.numeric(strs[length(strs)-1])){
		id==-1;
	}
	else{
		id= which(ids==as.numeric(strs[length(strs)-1]))
	}
}
	if(length(id)==1){
		ancestor[i] = id
	}
	else{ #look up ancestor of ancestor
		my=which(data[,2]==d[i,2])
		strs=strsplit(as.character(data[my,3]),'/')[[1]];
		id= which(as.numeric(data[,1])==as.numeric(strs[length(strs)-2]))
		ancestor[i]=which(d[,2]==data[id,2])
		print(ancestor[i])
#		ancestor[i]=-2;
	}
}








depth=matrix(0,n,1)
totalOffspring=matrix(0,n,1)#how much below
totalVol=matrix(0,n,1)#volume
for(j in 1:15){
	for(i in 1:n){
		if(ancestor[i]==-1){
			depth[i]=1
		}
		else{
			depth[i]=depth[ancestor[i]]+1
			
			
		}
		
		totalOffspring[i]=sum(ancestor==i)+sum(totalOffspring[ancestor==i])#how much below
		totalVol[i]=d$Volume..original.[i]+sum(totalVol[ancestor==i])#how much below
		if(sum(ancestor==i)==0&totalVol[i]<2){
			totalVol[i]=2;
		}
	}
	
}


#add dummies for each tip
for(i in 1:n){
	if(sum(ancestor==i)==0){
		##			add a dummy node, to preserve all info
			d=rbind(d,c(-1,-1,-1,'Dummy','DUM',1,1,0,NA))
			n=dim(d)[1]
			ancestor[n]=i;
			totalVol[n]=1;
			depth[n]=depth[i]+1
			d=rbind(d,c(-1,-1,-1,'Dummy','DUM',totalVol[i]-1,totalVol[i]-1,0,NA))
			n=dim(d)[1]
			ancestor[n]=i;
			depth[n]=depth[i]+1
			totalVol[n]=totalVol[i]-1;
	}
}




#get the difference between volums: height of the trees
diffVols=matrix(0,n,1)
for(i in 1:n){
	if(ancestor[i]!=-1){
		
	diffVols[i]=log10(totalVol[ancestor[i]])-log10(totalVol[i])
}
	else{
		diffVols[i]=log10(sum(totalVol[ancestor==-1]))-log10(totalVol[i])
	}
}



#check the diffVols add up to the same
#che=matrix(0,n,1)
#for(j in 1:15){
#	for(i in 1:n){
#		if(ancestor[i]==-1){
#			che[i]=diffVols[i]
#		}
#		else{
#			che[i]=che[ancestor[i]]+diffVols[i]
#		}
#		
#	}
#	
#}
#follow an edge up
#id=3
#s=0
#while(id%in%phyl$edge[,2]){
#	edgeid=which(phyl$edge[,2]==id)
#	print(phyl$edge.length[edgeid])
#	print(id)
#	id=phyl$edge[edgeid,1]
#	s=s+phyl$edge.length[edgeid]
#}
#s



#make a newick format

addPart <- function(id){
	if(sum(ancestor==id)==0){#leaf
#		return(strsplit(names[id],',')[[1]])
#		return(substring(names[id],1,3))
		return(paste(names[id],':',toString((diffVols[id])),sep=''))
	}
	string=c();
	for(i in which(ancestor==id))
		string=c(string,addPart(i))
	s=(string[1])
	if(length(string)>1)
	for(i in 2:(length(string)))
		s=paste(s,',',string[i],sep='')
	return(paste('(',s,')',(names[id]),':',toString((diffVols[id])),sep=''))
}

names=as.character(d$Abbreviation)
u=as.numeric(d$Parcellation)
#take out the commas
for(i in 1:length(names)){
	if(u[i]==1){
		names[i]=(strsplit(names[i],',')[[1]][1])
	}
	else
		names[i]=''
}

roots=which(ancestor==-1)
newick = addPart(roots[4])
#newick = addPart(31)
phyl=read.tree(text=paste(newick,';',sep=''))

#plot(phyl,type='radial',show.tip.label=T,cex=.6)




nodeUsed=''!=phyl$node.label#is this internal node in our parcellation?
ntip=length(phyl$tip.label)
isUsedEdge=nodeUsed[phyl$edge[,1]-ntip]

usedEdge=phyl$edge[,1]#belongs to which in our parcellation?
for(i in 1:length(usedEdge)){
	while(nodeUsed[usedEdge[i]-ntip]==0&&usedEdge[i]!=ntip+1){
		usedEdge[i]=phyl$edge[which(phyl$edge[,2]==usedEdge[i]),1]
	}
}
#MRCA in parcellation for each tip
usedNode=1:length(phyl$tip.label)
for(i in 1:length(usedNode)){
	usedNode[i]=phyl$edge[phyl$edge[,2]==i,1]
}
for(i in 1:length(usedNode)){
	while(nodeUsed[usedNode[i]-ntip]==0&&usedNode[i]!=ntip+1){
		usedNode[i]=phyl$edge[which(phyl$edge[,2]==usedNode[i]),1]
	}
}
#cheat, put node label on one of the tips
for(i in 1:phyl$Nnode){
	if(phyl$node.label[i]!=''){
#		find the middle descendant
		descendants=which(phyl$edge[,1]==i+ntip);
		id=phyl$edge[descendants[round(length(descendants)/2)],2]
		while(sum(phyl$edge[,1]==id)>0){
			descendants=which(phyl$edge[,1]==id);
			id=phyl$edge[descendants[round(length(descendants)/2)],2]
		}
		phyl$tip.label[id]=phyl$node.label[i]
	}
}
mypal = c("#1562A0", "#4ECDC4", "#1B676B", "#FF6B6B", "#C44D58")
library(RColorBrewer)
mypal=brewer.pal(11,'Set3')
mypaltips=brewer.pal(8,'Dark2')
mypal=c('#222222',rep(mypal,200))
mypaltips=c('#222222',rep(mypaltips,200))
mypaltips=mypal
#plot(phyl,type='fan',show.tip.label=T,tip.color=mypal[depth],edge.color=mypal[isUsedEdge*3+1],edge.width=2*isUsedEdge+1,label.offset=.3)
plot(phyl,show.tip.label=T,show.node.label=T,tip.color=mypal[usedNode-ntip],
		edge.color=mypal[usedEdge-ntip],label.offset=0,
		cex=.7,direction='downwards',srt=90,adj=.5)
plot(phyl,type='fan',show.tip.label=T,tip.color=mypaltips[usedNode-ntip],show.node.label=F,edge.color=mypal[usedEdge-ntip],edge.width=isUsedEdge*1+1,cex=.5,label.offset=.1,adj=.5,no.margin=T)

#plot(phyl,show.tip.label=T,tip.color=mypal[depth],edge.width=3*isUsedEdge+1)

#simple plot
#names[ancestor==237]=''
#names[237]=''
#plot(read.tree(text=paste(addPart(1187),';',sep='')),show.node.label=T)