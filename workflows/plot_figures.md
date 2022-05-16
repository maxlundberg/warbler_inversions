# Figure 1 

R code to generate Figure1:

```
fst_10kb=read.delim("freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.fst_10kb.out",header=T)
fst_10kb_maf=read.delim("freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.fst_10kb.maf0.1.out",header=T)

fst_variants=read.delim("freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.fst_div_region_adj_scaffolds.fst.2.out",header=T)


#Replace negative fst values with zero

replace_neg=function(x){
if(x<0){
	return(0)
	}
else{
	return(x)
	}
}


fst_variants[,3]<-sapply(fst_variants[,3],replace_neg)


#General plot parameters:
plot_low=-2.3
plot_high=1.2

col_scaff1="light blue"
col_scaff2="light green"
col_scaff3="lavender"


#Fst 10 kb plot parameters
fst_10kb_plot_col="blue"
fst_10kb_plot_lwd=2

#FST variant plot parameters
fst_variants_col="grey"
fst_variants_pch=16

#Specify distance between scaffolds
buffer=3E5

#Specify interval size to be plotted in adjacent scaffolds 
interval_size=2e6

#Southern ww
scaffold_top=-0.25
scaffold_bottom=-0.35
text_scaffold_pos=-0.4

#Northern ww
scaffold_north_top=-0.9
scaffold_north_bottom=-1

#Chiffchaff
scaffold_cc_top=-1.5
scaffold_cc_bottom=-1.6


#Outgroup 
scaffold_fl_top=-2.1
scaffold_fl_bottom=-2.2

#Synteny plot parameters
scaffold_synteny_top_ww_north=-0.55
scaffold_synteny_bottom_ww_north=-0.9

scaffold_synteny_top_cc=-1.2
scaffold_synteny_bottom_cc=-1.5

scaffold_synteny_top_fl=-1.8
scaffold_synteny_bottom_fl=-2.1

#synteny_col=rgb(255,165,0,max=255,alpha=125)
synteny_col="lightgoldenrod"
synteny_col_border="NA"
synteny_col_border=synteny_col

min_aln=2000

#################################################################################################################################

#Chr1

#################################################################################################################################



###########################################

#SOUTHERN WW

###########################################

#Add 2 Mb of neighbouring scaffolds + 2*buffer + length of Scaffold19,which is 11,852,059 bp
#That is: interval_size*2+buffer*2+11.852E6 = 16452000 bp

plot_start=-1.3E6
plot_end=16.452E6

##############

#Upstream scaffold - Scaffold12  14,402,903 bp - orientation: reversed
fst_10kb_maf_scaff=fst_10kb_maf[which(fst_10kb_maf$CHROM=="Scaffold12"),]

max_pos=2E6


#Plot FST for variants
fst_variants_scaff=fst_variants[which(fst_variants$CHROM=="Scaffold12"),]
fst_variants_scaff=fst_variants_scaff[which(fst_variants_scaff[,2]<=interval_size),]
plot(interval_size-fst_variants_scaff[,2],fst_variants_scaff[,3],xlim=c(plot_start,plot_end),ylim=c(plot_low,plot_high),axes=F,xlab="",ylab="",main="Chromosome 1",pch=fst_variants_pch,col=fst_variants_col)
axis(2,pos=-5E5,at=c(0:10)/10,cex.axis=0.7)
text(y=0.5,x=-1.9E6,"FST",srt=90,cex=0.8)

#Plot FST for 10 kb windows
fst_10kb_maf_scaff=fst_10kb_maf_scaff[which(fst_10kb_maf_scaff[,3]<=interval_size),]
points(interval_size-rowMeans(fst_10kb_maf_scaff[,2:3]),fst_10kb_maf_scaff[,5],type="l",col=fst_10kb_plot_col,lwd=fst_10kb_plot_lwd)

#Add scaffold segment
polygon(x=c(0,2E6,2E6,0),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col=col_scaff1,border=col_scaff1)
text(1E6,scaffold_bottom-0.1,"12 (-)",cex=0.7)

#Add tandem repeat array
#Scaffold12:13-66660 (1 based) -> add 1-66660
tr_start=interval_size-66660
tr_end=interval_size

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col="red",border="red")

##############

#Divergent region scaffold - Scaffold19 - length: 11,852,059 bp

cumulative_pos=buffer+interval_size

#To better show the two inverted intervals, we should plot this scaffold in the reverse orientation 


#Plot FST for variants 
fst_variants_scaff=fst_variants[which(fst_variants$CHROM=="Scaffold19"),]
points(11852059-fst_variants_scaff[,2]+cumulative_pos,fst_variants_scaff[,3],pch=fst_variants_pch,col=fst_variants_col)

#Plot FST for 10 kb windows
fst_10kb_maf_scaff=fst_10kb_maf[which(fst_10kb_maf$CHROM=="Scaffold19"),]
points(11852059-rowMeans(fst_10kb_maf_scaff[,2:3])+cumulative_pos,fst_10kb_maf_scaff[,5],type="l",col=fst_10kb_plot_col,lwd=fst_10kb_plot_lwd)


#Add scaffold segment
polygon(x=c(cumulative_pos,11852059+cumulative_pos,11852059+cumulative_pos,cumulative_pos),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col=col_scaff2,border=col_scaff2)
text(mean(c(1,11852059))+cumulative_pos,scaffold_bottom-0.1,"19 (-)",cex=0.7)


#Add tandem repeat intervals
#Start:1-48983 - 1 based
tr_start=11852059-48963+cumulative_pos
tr_end=11852059+cumulative_pos
polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col="red",border="red")
#End:11678266-11852059 ~ 174 kb - 1-based
tr_start=11852059-11852059+cumulative_pos
tr_end=11852059-11678265+cumulative_pos
polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col="red",border="red")
#Middle: 7978793-7980340 - 1-based (1548 bp)
tr_start=11852059-7980340+cumulative_pos
tr_end=11852059-7978793+cumulative_pos
polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col="red",border="red")

#Annotate the two inversion intervals 

#Interval nr 1: 58075-7968505 (based on mummer) - since the scaffold is reversed, in the plot, this will be the second interval
segments(11852059-7968505+cumulative_pos,1.2,11852059+cumulative_pos-58075,1.2)
text(mean(c(11852059-7968505,11852059))+cumulative_pos,1.28,"8.0 Mb",cex=0.7)
#Interval nr 2: 7978783-11844718 (based on mummer)
segments(11852059-7978783+cumulative_pos,1.1,11852059-11844718+cumulative_pos,1.1)
text(mean(c(11852059-7978783,11852059-11844718))+cumulative_pos,1.18,"3.9 Mb",cex=0.7)


##############

#Downstream scaffold (reverse orientation) - Scaffold11 - 47,735,347 bp

cumulative_pos=cumulative_pos+11852059+buffer

end_pos=47735347


#Need to fix boundary

#Plot FST for variants 
fst_variants_scaff=fst_variants[which(fst_variants$CHROM=="Scaffold11" & fst_variants$POS>=(end_pos-interval_size)),]
points(end_pos-fst_variants_scaff[,2]+cumulative_pos,fst_variants_scaff[,3],pch=fst_variants_pch,col=fst_variants_col)

#Plot FST for 10 kb windows
fst_10kb_maf_scaff=fst_10kb_maf[which(fst_10kb_maf$CHROM=="Scaffold11" & fst_10kb_maf$BIN_START>=(end_pos-interval_size)),]
fst_10kb_maf_scaff=fst_10kb_maf_scaff[which(fst_10kb_maf_scaff[,2]<=(end_pos+interval_size)),]
points(end_pos-rowMeans(fst_10kb_maf_scaff[,2:3])+cumulative_pos,fst_10kb_maf_scaff[,5],type="l",col=fst_10kb_plot_col,lwd=fst_10kb_plot_lwd)

#Add scaffold segment
start_scaffold=cumulative_pos
end_scaffold=cumulative_pos+interval_size
polygon(x=c(start_scaffold,end_scaffold,end_scaffold,start_scaffold),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col=col_scaff3,border=col_scaff3)
text(mean(c(start_scaffold,end_scaffold)),scaffold_bottom-0.1,"11 (-)",cex=0.7)

#TR:47,665,304-47,735,347 bp
tr_start=start_scaffold
tr_end=start_scaffold+end_pos-47665304
polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col="red",border="red")

##############

#Add sample label
text(y=mean(c(scaffold_bottom,scaffold_top)),x=-1E6,"southern",cex=0.8)


###########################################

#NORTHERN WW

###########################################


chr1_northern_mummer=read.delim("chr1.scaffolds.ww_north.vs.southern_ww.filt.coords.out",header=F)

#Filter based on alignment size
chr1_northern_mummer_filt=chr1_northern_mummer[which(chr1_northern_mummer[,5]>=min_aln & chr1_northern_mummer[,6]>=min_aln ),]


##############

#Upstream scaffold - Scaffold374 - 15,960,556 bp (reverse orientation relative outgroup)

#Scaffold12      68558   75397   Scaffold374     0       6843    0.995467        +


polygon(x=c(0,2E6,2E6,0),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col=col_scaff1,border=col_scaff1)
text(1E6,scaffold_north_bottom-0.1,"374 (-)",cex=0.7)

#Connect lines between this scaffold and the corresponding one in the southern genome


chr1_northern_mummer_filt_scaff=chr1_northern_mummer_filt[which(chr1_northern_mummer_filt[,10]=="Scaffold12" & chr1_northern_mummer_filt[,11]=="Scaffold374" & chr1_northern_mummer_filt[,4]<=interval_size),]


for(i in c(1:nrow(chr1_northern_mummer_filt_scaff))){

	input=chr1_northern_mummer_filt_scaff[i,]

	target_start=2E6-as.numeric(input[1])
	target_end=2E6-as.numeric(input[2])

	query_start=2E6-as.numeric(input[3])
	query_end=2E6-as.numeric(input[4])

	polygon(x=c(query_start,query_end,target_end,target_start),c(scaffold_synteny_bottom_ww_north,scaffold_synteny_bottom_ww_north,scaffold_synteny_top_ww_north,scaffold_synteny_top_ww_north),col=synteny_col,border=synteny_col_border)
}



#polygon(x=c(southern_start,southern_end,warbler_start,warbler_end),y=c(rep(southern_bottom,2),rep(southern_top,2)),col=warbler_color1,border=NA) 




##############

#Divergent region - Scaffold156 = 11,658,755 bp
cumulative_pos=2E6+buffer
end_scaffold=cumulative_pos+11658755

polygon(x=c(cumulative_pos,end_scaffold,end_scaffold,cumulative_pos),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col=col_scaff2,border=col_scaff2)
text(mean(c(cumulative_pos,end_scaffold)),scaffold_north_bottom-0.1,"156 (+)",cex=0.7)
 
#Scaffold19	58107	59385	Scaffold156	3743796	3745074	0.920188	+
#Scaffold19      7962885 7966543 Scaffold156     11653121        11656793        0.978403        +
#Scaffold19      7979025 7979994 Scaffold156     168     1141    0.821465        +
#Scaffold19      11683595        11683802        Scaffold156     3675467 3675674 0.710145        +

#TR at the start: Scaffold19:1-1952
tr_start=cumulative_pos+1
tr_end=cumulative_pos+1952

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col="red",border="red")

#TR in the middle - here we can include the gap (1-based)
#Scaffold156     3675482 3685135
#Scaffold156     3742954 3743449 

tr_start=cumulative_pos+3675482
tr_end=cumulative_pos+3743449

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col="red",border="red")

#Add gap in the middle: 3,685,149-3,742,948

gap_start=cumulative_pos+3685149
gap_end=cumulative_pos+3742948

polygon(x=c(gap_start,gap_end,gap_end,gap_start),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col="black",border="black")

#Connect lines between this scaffold and the corresponding one in the southern genome

chr1_northern_mummer_filt_scaff=chr1_northern_mummer_filt[which(chr1_northern_mummer_filt[,10]=="Scaffold19" & chr1_northern_mummer_filt[,11]=="Scaffold156"),]


for(i in c(1:nrow(chr1_northern_mummer_filt_scaff))){

	input=chr1_northern_mummer_filt_scaff[i,]

	target_start=11852059-as.numeric(input[1])+cumulative_pos
	target_end=11852059-as.numeric(input[2])+cumulative_pos

	query_start=as.numeric(input[3])+cumulative_pos
	query_end=as.numeric(input[4])+cumulative_pos

	polygon(x=c(query_start,query_end,target_end,target_start),c(scaffold_synteny_bottom_ww_north,scaffold_synteny_bottom_ww_north,scaffold_synteny_top_ww_north,scaffold_synteny_top_ww_north),col=synteny_col,border=synteny_col_border)
}




##############

#Downstream scaffold - Scaffold143 - 14,465,578 bp (positive orientation)
#Scaffold11      47661085        47665202        Scaffold143     0       4117    0.995142        -

scaffold_start=cumulative_pos+11658755+buffer
scaffold_end=scaffold_start+2E6

polygon(x=c(scaffold_start,scaffold_end,scaffold_end,scaffold_start),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col=col_scaff3,border=col_scaff3)
text(mean(c(scaffold_start,scaffold_end)),scaffold_north_bottom-0.1,"143 (+)",cex=0.7)


chr1_northern_mummer_filt_scaff=chr1_northern_mummer_filt[which(chr1_northern_mummer_filt[,10]=="Scaffold11" & chr1_northern_mummer_filt[,11]=="Scaffold143" & chr1_northern_mummer_filt[,3]<=interval_size),]


for(i in c(1:nrow(chr1_northern_mummer_filt_scaff))){

	input=chr1_northern_mummer_filt_scaff[i,]

	target_start=47735347-as.numeric(input[1])+interval_size+11852059+2*buffer
	target_end=47735347-as.numeric(input[2])+interval_size+11852059+2*buffer

	query_start=as.numeric(input[3])+scaffold_start
	query_end=as.numeric(input[4])+scaffold_start

	polygon(x=c(query_start,query_end,target_end,target_start),c(scaffold_synteny_bottom_ww_north,scaffold_synteny_bottom_ww_north,scaffold_synteny_top_ww_north,scaffold_synteny_top_ww_north),col=synteny_col,border=synteny_col_border)
}


##############

#Add sample label

text(y=mean(c(scaffold_north_bottom,scaffold_north_top)),x=-1E6,"northern",cex=0.8)


###########################################

#CHIFFCHAFF

###########################################

chr1_cc_mummer=read.delim("chr1.scaffolds.cc.vs.southern_ww.filt.coords.out",header=F)

#Filter based on alignment size

chr1_cc_mummer_filt=chr1_cc_mummer[which(chr1_cc_mummer[,5]>=min_aln & chr1_cc_mummer[,6]>=min_aln),]


##############

#Upstream scaffold - reverse orientation - ptg000006l 62,520,444 bp

polygon(x=c(0,2E6,2E6,0),y=c(scaffold_cc_bottom,scaffold_cc_bottom,scaffold_cc_top,scaffold_cc_top),col=col_scaff1,border=col_scaff1)
text(1E6,scaffold_cc_bottom-0.1,"6l (-)",cex=0.7)

#Add tandem repeat (1-based pos)
#ptg000006l 1 613112 613112
#ptg000006l 613098 613151 54   


tr_start=2E6
tr_end=2E6-613151

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_cc_bottom,scaffold_cc_bottom,scaffold_cc_top,scaffold_cc_top),col="red",border="red")


#Plot synteny

chr1_cc_mummer_filt_scaff=chr1_cc_mummer_filt[which(chr1_cc_mummer_filt[,10]=="Scaffold12" & chr1_cc_mummer_filt[,11]=="ptg000006l" & chr1_cc_mummer_filt[,4]<=interval_size),]


for(i in c(1:nrow(chr1_cc_mummer_filt_scaff))){

	input=chr1_cc_mummer_filt_scaff[i,]

	target_start=2E6-as.numeric(input[1])
	target_end=2E6-as.numeric(input[2])

	query_start=2E6-as.numeric(input[3])
	query_end=2E6-as.numeric(input[4])

	polygon(x=c(query_start,query_end,target_end,target_start),c(scaffold_synteny_bottom_cc,scaffold_synteny_bottom_cc,scaffold_synteny_top_cc,scaffold_synteny_top_cc),col=synteny_col,border=synteny_col_border)
}




##############

#Divergent region scaffold - reverse orientation ptg000007l - 19,027,237 bp

cumulative_pos=2E6+buffer

#polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_cc_bottom,scaffold_cc_bottom,scaffold_cc_top,scaffold_cc_top),col="red",border=NA)


#Plot scaffold segment 
polygon(x=c(cumulative_pos,cumulative_pos+14E6,cumulative_pos+14E6,cumulative_pos),y=c(scaffold_cc_bottom,scaffold_cc_bottom,scaffold_cc_top,scaffold_cc_top),col=col_scaff2,border=col_scaff2)
text(mean(c(cumulative_pos,cumulative_pos+14E6)),scaffold_cc_bottom-0.1,"7l (-)",cex=0.7)

#Add TR at the end of the div region: 6512119-6952665 (1-based) ~440 kb

tr_start=cumulative_pos+19027237-6512120
tr_end=cumulative_pos+19027237-6952665

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_cc_bottom,scaffold_cc_bottom,scaffold_cc_top,scaffold_cc_top),col="red",border="red")

#Add TR in the central part of the div region: 14900329-15170453 (1-based) ~ 270 kb

tr_start=cumulative_pos+19027237-14900328
tr_end=cumulative_pos+19027237-15170453

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_cc_bottom,scaffold_cc_bottom,scaffold_cc_top,scaffold_cc_top),col="red",border="red")


#Add TR at the start of the div region. Here there are two close intervals: 
#1) 18874730 - 18983007 (1 based) ~ 108 kb

tr_start=cumulative_pos+19027237-18874730
tr_end=cumulative_pos+19027237-18983007

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_cc_bottom,scaffold_cc_bottom,scaffold_cc_top,scaffold_cc_top),col="red",border=NA)

#2 18989765 - 19027234 (1 based) - 37.5 kb 

tr_start=cumulative_pos+19027237-18989765
tr_end=cumulative_pos+19027237-19027234

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_cc_bottom,scaffold_cc_bottom,scaffold_cc_top,scaffold_cc_top),col="red",border=NA)



#Plot synteny 

#first part to the divergent region

chr1_cc_mummer_filt_scaff=chr1_cc_mummer_filt[which(chr1_cc_mummer_filt[,10]=="Scaffold19" & chr1_cc_mummer_filt[,11]=="ptg000007l"),]


for(i in c(1:nrow(chr1_cc_mummer_filt_scaff))){

	input=chr1_cc_mummer_filt_scaff[i,]

	target_start=cumulative_pos+11852059-as.numeric(input[1])
	target_end=cumulative_pos+11852059-as.numeric(input[2])

	query_start=cumulative_pos+19027237-as.numeric(input[3])
	query_end=cumulative_pos+19027237-as.numeric(input[4])

	polygon(x=c(query_start,query_end,target_end,target_start),c(scaffold_synteny_bottom_cc,scaffold_synteny_bottom_cc,scaffold_synteny_top_cc,scaffold_synteny_top_cc),col=synteny_col,border=synteny_col_border)
}




#second part to the downstream scaffold

#We should start after the end of the tandem repeat array in the chiffchaff and continue to interval_size - tandem repeat array length

chr1_cc_mummer_filt_scaff=chr1_cc_mummer_filt[which(chr1_cc_mummer_filt[,10]=="Scaffold11" & chr1_cc_mummer_filt[,11]=="ptg000007l" & chr1_cc_mummer_filt[,1]>=46E6),]


for(i in c(1:nrow(chr1_cc_mummer_filt_scaff))){

	input=chr1_cc_mummer_filt_scaff[i,]

	target_start=47735347-as.numeric(input[1])+interval_size+11852059+2*buffer
	target_end=47735347-as.numeric(input[2])+interval_size+11852059+2*buffer

	query_start=cumulative_pos+19027237-as.numeric(input[3])
	query_end=cumulative_pos+19027237-as.numeric(input[4])

	polygon(x=c(query_start,query_end,target_end,target_start),c(scaffold_synteny_bottom_cc,scaffold_synteny_bottom_cc,scaffold_synteny_top_cc,scaffold_synteny_top_cc),col=synteny_col,border=synteny_col_border)
}


##############

#Add sample label
text(y=mean(c(scaffold_cc_bottom,scaffold_cc_top)),x=-1E6,"chiffchaff",cex=0.8)



###########################################

#FLYCATCHER

###########################################

#Add chr1 flycatcher - orientation will be positive

#end_scaffold=plot_end

polygon(x=c(0,plot_end,plot_end,0),y=c(scaffold_fl_bottom,scaffold_fl_bottom,scaffold_fl_top,scaffold_fl_top),col=col_scaff2,border=NA)
text(plot_end/2,scaffold_fl_bottom-0.1,"Chromosome 1:62.5-78.6 Mb (+)",cex=0.7)

#Synteny

chr1_flycatcher_mummer=read.delim("chr1.scaffolds.ww_south.vs.flycatcher.delta.filt.coords.out",header=F)

chr1_flycatcher_mummer_filt=chr1_flycatcher_mummer[which(chr1_flycatcher_mummer[,5]>=min_aln & chr1_flycatcher_mummer[,6]>=min_aln),]

#Upstream scaffold
chr1_flycatcher_mummer_filt_scaff=chr1_flycatcher_mummer_filt[which(chr1_flycatcher_mummer_filt[,11]=="Scaffold12" & chr1_flycatcher_mummer_filt[,4]<=2E6),]


#62501826        62502474        2000227 1999577 649     651     91.44   1       -1      1       Scaffold12



for(i in c(1:nrow(chr1_flycatcher_mummer_filt_scaff))){

	input=chr1_flycatcher_mummer_filt_scaff[i,]

	target_start=2E6-as.numeric(input[3])
	target_end=2E6-as.numeric(input[4])

	query_start=as.numeric(input[1])-62.50E6
	query_end=as.numeric(input[2])-62.50E6

	polygon(x=c(query_start,query_end,target_end,target_start),c(scaffold_synteny_bottom_fl,scaffold_synteny_bottom_fl,scaffold_synteny_top_fl,scaffold_synteny_top_fl),col=synteny_col,border=synteny_col_border)
}

#Divergent scaffold

chr1_flycatcher_mummer_filt_scaff=chr1_flycatcher_mummer_filt[which(chr1_flycatcher_mummer_filt[,11]=="Scaffold19"),]


#62501826        62502474        2000227 1999577 649     651     91.44   1       -1      1       Scaffold12



for(i in c(1:nrow(chr1_flycatcher_mummer_filt_scaff))){

	input=chr1_flycatcher_mummer_filt_scaff[i,]

	target_start=2E6+buffer+11852059-as.numeric(input[3])
	target_end=2E6+buffer+11852059-as.numeric(input[4])

	query_start=as.numeric(input[1])-62.50E6
	query_end=as.numeric(input[2])-62.50E6

	polygon(x=c(query_start,query_end,target_end,target_start),c(scaffold_synteny_bottom_fl,scaffold_synteny_bottom_fl,scaffold_synteny_top_fl,scaffold_synteny_top_fl),col=synteny_col,border=synteny_col_border)
}


#Downstream scaffold

chr1_flycatcher_mummer_filt_scaff=chr1_flycatcher_mummer_filt[which(chr1_flycatcher_mummer_filt[,11]=="Scaffold11" & chr1_flycatcher_mummer_filt[,4]>=45.73E6), ,]

for(i in c(1:nrow(chr1_flycatcher_mummer_filt_scaff))){

	input=chr1_flycatcher_mummer_filt_scaff[i,]

	target_start=47735347-as.numeric(input[3])+interval_size+11852059+2*buffer
	target_end=47735347-as.numeric(input[4])+interval_size+11852059+2*buffer

	query_start=as.numeric(input[1])-62.50E6
	query_end=as.numeric(input[2])-62.50E6

	polygon(x=c(query_start,query_end,target_end,target_start),c(scaffold_synteny_bottom_fl,scaffold_synteny_bottom_fl,scaffold_synteny_top_fl,scaffold_synteny_top_fl),col=synteny_col,border=synteny_col_border)
}

#Add sample label
text(y=mean(c(scaffold_fl_bottom,scaffold_fl_top)),x=-1E6,"flycatcher",cex=0.8)




#################################################################################################################################

#Chr3

#################################################################################################################################

plot_start=-1.3E6
plot_end=17.6E6

###########################################

#SOUTHERN WW

###########################################


##############
#Upstream scaffold (positive orientation) - Scaffold38 - 35,991,252 bp

cumulative_pos=end_scaffold+buffer


fst_10kb_maf_scaff=fst_10kb_maf[which(fst_10kb_maf$CHROM=="Scaffold38" & fst_10kb_maf$BIN_START>=(35991252-interval_size)),]
fst_variants_scaff=fst_variants[which(fst_variants$CHROM=="Scaffold38" & fst_variants$POS>=(35991252-interval_size)),]

#scaffold_end=35991252 

#Plot FST for variants
plot(as.numeric(paste(fst_variants_scaff[,2]))-(35991252-interval_size),as.numeric(paste(fst_variants_scaff[,3])),xlim=c(plot_start,plot_end),ylim=c(plot_low,plot_high),axes=F,xlab="",ylab="",main="Chromosome 3",pch=fst_variants_pch,col=fst_variants_col)

#Plot FST for 10 kb windows
points(rowMeans(fst_10kb_maf_scaff[,2:3])-(35991252-interval_size),fst_10kb_maf_scaff[,5],type="l",col=fst_10kb_plot_col,lwd=fst_10kb_plot_lwd)

axis(2,pos=-5E5,at=c(0:10)/10,cex.axis=0.7)

#text(x=6E6,y=-0.6,"position (Mb)")
text(y=0.5,x=-1.9E6,"FST",srt=90,cex=0.8)


#Plot scaffold

start_scaffold=0
end_scaffold=interval_size
polygon(x=c(start_scaffold,end_scaffold,end_scaffold,start_scaffold),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col=col_scaff1,border=col_scaff1)
text(1E6,scaffold_bottom-0.1,"38 (+)",cex=0.7)

#plot_end+cumulative_pos

#Plot TR - 35897021-35991244 (1-based) ~ 94224

tr_start=35897021-(35991252-interval_size)
tr_end=35991244-(35991252-interval_size)

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col="red",border="red")

text(y=mean(c(scaffold_bottom,scaffold_top)),x=-1E6,"southern",cex=0.8)


##############
#Divergent scaffold - Scaffold61 - orientation: reversed - 69,302,477 bp

cumulative_pos=interval_size+buffer

fst_10kb_maf_scaff=fst_10kb_maf[which(fst_10kb_maf$CHROM=="Scaffold61"),]

scaffold_start=53.89E6

#Plot FST for variants
fst_variants_scaff=fst_variants[which(fst_variants$CHROM=="Scaffold61"),]
fst_variants_scaff=fst_variants_scaff[which(fst_variants_scaff[,2]>=scaffold_start),]
points(69302477-fst_variants_scaff[,2]+cumulative_pos,fst_variants_scaff[,3],col=fst_variants_col,pch=fst_variants_pch)
#axis(1,pos=-0.2,at=c(0:15)*1E6,labels=c(55:70),cex=0.8)

#Plot FST for 10 kb windows
fst_10kb_maf_scaff=fst_10kb_maf_scaff[which(fst_10kb_maf_scaff[,3]>=scaffold_start),]
points(69302477-rowMeans(fst_10kb_maf_scaff[,2:3])+cumulative_pos,fst_10kb_maf_scaff[,5],type="l",col=fst_10kb_plot_col,lwd=fst_10kb_plot_lwd)


#Add scaffold segment
start_scaffold=interval_size+buffer
end_scaffold=start_scaffold+69302477-53.89E6

polygon(x=c(start_scaffold,end_scaffold,end_scaffold,start_scaffold),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col=col_scaff2,border=col_scaff2)
text(mean(c(start_scaffold,end_scaffold)),scaffold_bottom-0.1,"61 (-)",cex=0.7)


#Add TR

#There are two intervals at the end - 
#1) Scaffold61:69290087-69302465 (1-based) ~ 12379 bp
tr_start=69302477-69290087+cumulative_pos
tr_end=69302477-69302465+cumulative_pos
polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col="red",border="red")

#2) Scaffold61:69204988-69283336 (1-based) ~ 78349 bp
tr_start=69302477-69204988+cumulative_pos
tr_end=69302477-69283336+cumulative_pos

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col="red",border="red")


#TR at the other end of the divergent region: 55894452-56078109 (1-based) ~ 183,658 bp
tr_start=69302477-55894451+cumulative_pos
tr_end=69302477-56078109+cumulative_pos

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col="red",border="red")


#Add line to show the divergent region - we could define it based on the synteny
#56069986-69205658
segments(69302477-69205658+cumulative_pos,1.1,69302477-56069986+cumulative_pos,1.1)
text(mean(c(69302477-69205658,69302477-56069986))+cumulative_pos,1.18,"13.1 Mb",cex=0.7)


###########################################

#NORTHERN WW

###########################################

scaffold_synteny_top=-0.4
scaffold_synteny_bottom=-0.65
#synteny_col=rgb(255,165,0,max=255,alpha=125)

##############

#Upstream scaffold - Scaffold139 (reverse) - 32,231,478 bp

start_scaffold=0
end_scaffold=interval_size

polygon(x=c(start_scaffold,end_scaffold,end_scaffold,start_scaffold),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col=col_scaff1,border=col_scaff1)
text(1E6,scaffold_north_bottom-0.1,"139 (-)",cex=0.7)


#Add TR - 1-3135 (1-based)

tr_start=interval_size-3135
tr_end=interval_size

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col="red",border="red")

#Synteny 
chr3_northern_mummer=read.delim("chr3.scaffolds.ww_north.vs.southern_ww.filt.coords.out",header=F)
chr3_northern_mummer_filt=chr3_northern_mummer[which(chr3_northern_mummer[,5]>=min_aln & chr3_northern_mummer[,6]>=min_aln),]

#The [,4] match looks better, because there is a large synteny block at the end
chr3_northern_mummer_filt_scaff=chr3_northern_mummer_filt[which(chr3_northern_mummer_filt[,10]=="Scaffold38" & chr3_northern_mummer_filt[,11]=="Scaffold139" & chr3_northern_mummer_filt[,4]<=interval_size),]


for(i in c(1:nrow(chr3_northern_mummer_filt_scaff))){

	input=chr3_northern_mummer_filt_scaff[i,]

	target_start=as.numeric(input[1])-(35991252-interval_size)
	target_end=as.numeric(input[2])-(35991252-interval_size)

	query_start=interval_size-as.numeric(input[3])
	query_end=interval_size-as.numeric(input[4])

	polygon(x=c(query_start,query_end,target_end,target_start),c(scaffold_synteny_bottom_ww_north,scaffold_synteny_bottom_ww_north,scaffold_synteny_top_ww_north,scaffold_synteny_top_ww_north),col=synteny_col,border=synteny_col_border)
}

text(y=mean(c(scaffold_north_bottom,scaffold_north_top)),x=-1E6,"northern",cex=0.8)

##############

#Divergent scaffold - Scaffold29b - 13,167,625 bp
cumulative_pos=2E6+buffer

start_scaffold=cumulative_pos
end_scaffold=cumulative_pos+13167625

polygon(x=c(start_scaffold,end_scaffold,end_scaffold,start_scaffold),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col=col_scaff2,border=col_scaff2)
text(mean(c(start_scaffold,end_scaffold)),scaffold_north_bottom-0.1,"29b (-)",cex=0.7)

#Add TR at the start of the scaffold (end of region) - 55759505-55767600 (1-based) ~8 kb

tr_start=13167625-(55759505-55759505)+cumulative_pos
tr_end=13167625-(55767600-55759505)+cumulative_pos

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col="red",border="red")



#Synteny
chr3_northern_mummer_filt_scaff=chr3_northern_mummer_filt[which(chr3_northern_mummer_filt[,10]=="Scaffold61" & chr3_northern_mummer_filt[,11]=="Scaffold29b"),]


for(i in c(1:nrow(chr3_northern_mummer_filt_scaff))){

	input=chr3_northern_mummer_filt_scaff[i,]

	target_start=69302477-as.numeric(input[1])+cumulative_pos
	target_end=69302477-as.numeric(input[2])+cumulative_pos

	query_start=13167625-as.numeric(input[3])+cumulative_pos
	query_end=13167625-as.numeric(input[4])+cumulative_pos

	polygon(x=c(query_start,query_end,target_end,target_start),c(scaffold_synteny_bottom_ww_north,scaffold_synteny_bottom_ww_north,scaffold_synteny_top_ww_north,scaffold_synteny_top_ww_north),col=synteny_col,border=synteny_col_border)
}


##############

#Downstream scaffold - Scaffold29 (reverse orientation) - 55,738,393 bp

cumulative_pos=cumulative_pos+13167625+buffer

polygon(x=c(cumulative_pos,cumulative_pos+interval_size,cumulative_pos+interval_size,cumulative_pos),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col=col_scaff3,border=col_scaff3)
text(mean(c(cumulative_pos,cumulative_pos+interval_size)),scaffold_north_bottom-0.1,"29 (-)",cex=0.7)

#Add TR - 55733727-55738385 (1-based) ~ 4.7 kb

tr_start=55738393-55733727+cumulative_pos 
tr_end=55738393-55738385+cumulative_pos

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col="red",border="red")


#Synteny 

#Connect lines between this scaffold and the corresponding one in the southern genome

chr3_northern_mummer_filt_scaff=chr3_northern_mummer_filt[which(chr3_northern_mummer_filt[,10]=="Scaffold61" & chr3_northern_mummer_filt[,11]=="Scaffold29" & chr3_northern_mummer_filt[,3]>=(55738393-interval_size)),]


for(i in c(1:nrow(chr3_northern_mummer_filt_scaff))){

	input=chr3_northern_mummer_filt_scaff[i,]

	target_start=69302477-as.numeric(input[1])+buffer+interval_size
	target_end=69302477-as.numeric(input[2])+buffer+interval_size

	query_start=55738393-as.numeric(input[3])+cumulative_pos
	query_end=55738393-as.numeric(input[4])+cumulative_pos

	polygon(x=c(query_start,query_end,target_end,target_start),c(scaffold_synteny_bottom_ww_north,scaffold_synteny_bottom_ww_north,scaffold_synteny_top_ww_north,scaffold_synteny_top_ww_north),col=synteny_col,border=synteny_col_border)
}






###########################################

#Chiffchaff

###########################################



chr3_cc_mummer=read.delim("chr3.scaffolds.cc.vs.southern_ww.filt.coords.out",header=F)

#Filter based on alignment size
chr3_cc_mummer_filt=chr3_cc_mummer[which(chr3_cc_mummer[,5]>=min_aln & chr3_cc_mummer[,5]>=min_aln),]

#################

#Divergent region/downstream scaffold ptg000040l - 29,328,574 bp -reverse orientation

start=0
end=15.5E6

polygon(x=c(start,end,end,start),y=c(scaffold_cc_bottom,scaffold_cc_bottom,scaffold_cc_top,scaffold_cc_top),col=col_scaff2,border=col_scaff2)
text(end/2,scaffold_cc_bottom-0.1,"40l (+)",cex=0.7)

#Add TR

#ptg000040l:15414192-15858620 (1-based) - 444429 bp


tr_start=15414192-(29328574-15.5E6)
tr_end=15858620-(29328574-15.5E6)

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_cc_bottom,scaffold_cc_bottom,scaffold_cc_top,scaffold_cc_top),col="red",border="red")


#ptg000040l:29049553-29328563 (1-based) - 279,011 bp

tr_start=29049552-(29328574-15.5E6)
tr_end=29328563-(29328574-15.5E6)

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_cc_bottom,scaffold_cc_bottom,scaffold_cc_top,scaffold_cc_top),col="red",border="red")

#Final TR at the very end: 29269249-29328563 (1-based) - 59,315 bp

tr_start=29269249-(29328574-15.5E6)
tr_end=29328563-(29328574-15.5E6)

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_cc_bottom,scaffold_cc_bottom,scaffold_cc_top,scaffold_cc_top),col="red",border="red")



#Synteny

#To scaffold38 in the southern willow warbler assembly

chr3_cc_mummer_filt_scaff=chr3_cc_mummer_filt[which(chr3_cc_mummer_filt[,10]=="Scaffold38" & chr3_cc_mummer_filt[,11]=="ptg000040l" & chr3_cc_mummer_filt[,3]>=(29328574-15.5E6) & chr3_cc_mummer_filt[,3]<(29328574-15.5E6+interval_size) & chr3_cc_mummer_filt[,1]>=(35991252-interval_size)),]

for(i in c(1:nrow(chr3_cc_mummer_filt_scaff))){

	input=chr3_cc_mummer_filt_scaff[i,]

	target_start=as.numeric(input[1])-(35991252-interval_size)
	target_end=as.numeric(input[2])-(35991252-interval_size)

	query_start=as.numeric(input[3])-(29328574-15.5E6)
	query_end=as.numeric(input[4])-(29328574-15.5E6)

	polygon(x=c(query_start,query_end,target_end,target_start),c(scaffold_synteny_bottom_cc,scaffold_synteny_bottom_cc,scaffold_synteny_top_cc,scaffold_synteny_top_cc),col=synteny_col,,border=synteny_col_border)
}


#To div region on Scaffold61 in the southern willow warbler assembly

chr3_cc_mummer_filt_scaff=chr3_cc_mummer_filt[which(chr3_cc_mummer_filt[,10]=="Scaffold61" & chr3_cc_mummer_filt[,11]=="ptg000040l" & chr3_cc_mummer_filt[,1]>56E6),]

for(i in c(1:nrow(chr3_cc_mummer_filt_scaff))){

	input=chr3_cc_mummer_filt_scaff[i,]

	target_start=69302477-as.numeric(input[1])+interval_size+buffer
	target_end=69302477-as.numeric(input[2])+interval_size+buffer

	query_start=as.numeric(input[3])-(29328574-15.5E6)
	query_end=as.numeric(input[4])-(29328574-15.5E6)

	polygon(x=c(query_start,query_end,target_end,target_start),c(scaffold_synteny_bottom_cc,scaffold_synteny_bottom_cc,scaffold_synteny_top_cc,scaffold_synteny_top_cc),col=synteny_col,border=synteny_col_border)
}


text(y=mean(c(scaffold_cc_bottom,scaffold_cc_top)),x=-1E6,"chiffchaff",cex=0.8)

#################

#Downstream scaffold - ptg000026l - 56,453,846 - positive orientation

start=15.5E6+buffer
end=start+interval_size

polygon(x=c(start,end,end,start),y=c(scaffold_cc_bottom,scaffold_cc_bottom,scaffold_cc_top,scaffold_cc_top),col=col_scaff3,border=col_scaff3)
text(mean(c(start,end)),scaffold_cc_bottom-0.1,"26l (+)",cex=0.7)

#Add TR - 1-251924 (1-based) ~ 252 kb

tr_start=start
tr_end=start+251924

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_cc_bottom,scaffold_cc_bottom,scaffold_cc_top,scaffold_cc_top),col="red",border="red")


#Synteny

chr3_cc_mummer_filt_scaff=chr3_cc_mummer_filt[which(chr3_cc_mummer_filt[,10]=="Scaffold61" & chr3_cc_mummer_filt[,11]=="ptg000026l" & chr3_cc_mummer_filt[,4]<=interval_size),]



for(i in c(1:nrow(chr3_cc_mummer_filt_scaff))){

	input=chr3_cc_mummer_filt_scaff[i,]

	target_start=69302477-as.numeric(input[1])+interval_size+buffer
	target_end=69302477-as.numeric(input[2])+interval_size+buffer

	query_start=as.numeric(input[3])+start
	query_end=as.numeric(input[4])+start

	polygon(x=c(query_start,query_end,target_end,target_start),c(scaffold_synteny_bottom_cc,scaffold_synteny_bottom_cc,scaffold_synteny_top_cc,scaffold_synteny_top_cc),col=synteny_col,border=synteny_col_border)
}





###########################################

#Flycatcher

###########################################

#Add chr3 flycatcher - orientation will be positive

#end_scaffold=plot_end

polygon(x=c(0,plot_end,plot_end,0),y=c(scaffold_fl_bottom,scaffold_fl_bottom,scaffold_fl_top,scaffold_fl_top),col=col_scaff2,border=col_scaff2)
text(plot_end/2,scaffold_fl_bottom-0.1,"Chromosome 3: 42.4-60.0 Mb (+)",cex=0.7)

#Synteny

chr3_flycatcher_mummer=read.delim("chr3.scaffolds.ww_south.vs.flycatcher.delta.filt.coords.out",header=F)
min_aln=2000

chr3_flycatcher_mummer_filt=chr3_flycatcher_mummer[which(chr3_flycatcher_mummer[,5]>=min_aln & chr3_flycatcher_mummer[,6]>=min_aln),]

scaffold_synteny_top=-1.25
scaffold_synteny_bottom=-1.45


#Upstream scaffold

Scaffold38_size=35991252

#tr_start=35897020-(35991252-interval_size)
#tr_end=35991244-(35991252-interval_size)



chr3_flycatcher_mummer_filt_scaff=chr3_flycatcher_mummer_filt[which(chr3_flycatcher_mummer_filt[,11]=="Scaffold38" & chr3_flycatcher_mummer_filt[,4]>=(Scaffold38_size-interval_size)),]


#62501826        62502474        2000227 1999577 649     651     91.44   1       -1      1       Scaffold12



for(i in c(1:nrow(chr3_flycatcher_mummer_filt_scaff))){

	input=chr3_flycatcher_mummer_filt_scaff[i,]

	target_start=as.numeric(input[3])-(Scaffold38_size-interval_size)
	target_end=as.numeric(input[4])-(Scaffold38_size-interval_size)

	query_start=as.numeric(input[1])-42.38E6
	query_end=as.numeric(input[2])-42.38E6

	polygon(x=c(query_start,query_end,target_end,target_start),c(scaffold_synteny_bottom_fl,scaffold_synteny_bottom_fl,scaffold_synteny_top_fl,scaffold_synteny_top_fl),col=synteny_col,border=synteny_col_border)
}

#Divergent scaffold

chr3_flycatcher_mummer_filt_scaff=chr3_flycatcher_mummer_filt[which(chr3_flycatcher_mummer_filt[,11]=="Scaffold61" & chr3_flycatcher_mummer_filt[,2]<=plot_end+42.38E6),]

for(i in c(1:nrow(chr3_flycatcher_mummer_filt_scaff))){

	input=chr3_flycatcher_mummer_filt_scaff[i,]


	target_start=69302477-as.numeric(input[3])+interval_size+buffer
	target_end=69302477-as.numeric(input[4])+interval_size+buffer

	query_start=as.numeric(input[1])-42.38E6
	query_end=as.numeric(input[2])-42.38E6

	polygon(x=c(query_start,query_end,target_end,target_start),c(scaffold_synteny_bottom_fl,scaffold_synteny_bottom_fl,scaffold_synteny_top_fl,scaffold_synteny_top_fl),col=synteny_col,border=synteny_col_border)
}



text(y=mean(c(scaffold_fl_bottom,scaffold_fl_top)),x=-1E6,"flycatcher",cex=0.8)




#################################################################################################################################

#Chr5

#################################################################################################################################

plot_start=-0.65E6
plot_end=5E6


fst_10kb_maf_scaff=fst_10kb_maf[which(fst_10kb_maf$CHROM=="Scaffold0" & fst_10kb_maf$BIN_END<=plot_end),]
fst_variants_scaff=fst_variants[which(fst_variants$CHROM=="Scaffold0" & as.numeric(paste(fst_variants$POS))<=plot_end),]


#Plot FST for variants 

plot(as.numeric(paste(fst_variants_scaff[,2])),as.numeric(paste(fst_variants_scaff[,3])),xlim=c(plot_start,plot_end),ylim=c(plot_low,plot_high),axes=F,xlab="",ylab="",main="Chromosome 5",pch=fst_variants_pch,col=fst_variants_col)
#axis(1,pos=-0.2,at=c(0:15)*1E6,labels=c(55:70),cex=0.8)
axis(2,pos=-2E5,at=c(0:10)/10,cex.axis=0.7)

#text(x=6E6,y=-0.6,"position (Mb)")
text(y=0.5,x=-0.8E6,"FST",srt=90,cex=0.8)

#Plot fst for 10 kb windows
#fst_10kb_maf_scaff=fst_10kb_maf_scaff[which(fst_10kb_maf_scaff[,3]>=55E6),]
points(rowMeans(fst_10kb_maf_scaff[,2:3]),fst_10kb_maf_scaff[,5],type="l",col=fst_10kb_plot_col,lwd=fst_10kb_plot_lwd)


#Plot scaffold segment 

polygon(x=c(0,plot_end,plot_end,0),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col=col_scaff2,border=col_scaff2)
text(mean(c(plot_start,plot_end)),scaffold_bottom-0.1,"0 (+)",cex=0.7)


#TR at the start 212054-328479 (1-based) - 116 kb

tr_start=212054
tr_end=328479

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col="red",border="red")


#Add lines representing the two plateaus

#First plateau -starts with the first high (>=0.7) FST variant at 207,169 bp and continues to about 3915477 bp based on synteny

segments(207169,1.1,3915477,1.1)
text(mean(c(207169,3915477)),1.18,"3.7 Mb",cex=0.7)

#Second plateau - starts around 3917392 and ends around   4361809
segments(3917392,1.2,4361809,1.2)
text(mean(c(3917392,4361809)),1.28,"0.4 Mb",cex=0.7)


#Synteny in the southern genom
#4358697  4361809

#Add sample label
text(y=mean(c(scaffold_bottom,scaffold_top)),x=-0.5E6,"southern",cex=0.8)

###########################################

#NORTHERN WW

###########################################


#Scaffold68 - reverse orientation  -  4,634,130 bp

polygon(x=c(0,4634130,4634130,0),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col=col_scaff2,border=col_scaff2)
text(4634130/2,scaffold_north_bottom-0.1,"68 (-)",cex=0.7) 
#Synteny

chr5_northern_mummer=read.delim("chr5.scaffolds.ww_north.vs.southern_ww.filt.coords.out",header=F)

#Filter based on alignment size. There are some short intervals mapping further away in the southern genome. 

chr5_northern_mummer_filt=chr5_northern_mummer[which(chr5_northern_mummer[,5]>=min_aln & chr5_northern_mummer[,6]>=min_aln & chr5_northern_mummer[,1]<=5E6),]

#synteny_col=rgb(255,165,0,max=255,alpha=125)


for(i in c(1:nrow(chr5_northern_mummer_filt))){

	input=chr5_northern_mummer_filt[i,]

	target_start=as.numeric(input[1])
	target_end=as.numeric(input[2])

	query_start=4634130-as.numeric(input[3])
	query_end=4634130-as.numeric(input[4])

	synteny_col_chr5=synteny_col
	synteny_col_border_chr5=synteny_col
	
	if(target_start>=3917392 & target_end<=4361809){
		#synteny_col_chr5=rgb(255,0,0,max=255,alpha=125)
		synteny_col_chr5="coral"
		synteny_col_border_chr5="coral"
		}


	
	polygon(x=c(query_start,query_end,target_end,target_start),c(scaffold_synteny_bottom_ww_north,scaffold_synteny_bottom_ww_north,scaffold_synteny_top_ww_north,scaffold_synteny_top_ww_north),col=synteny_col_chr5,border=synteny_col_border_chr5)
}


#Add gap

gap_start=4634130-4046153
gap_end=4634130-4087200

polygon(x=c(gap_start,gap_end,gap_end,gap_start),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col="black",border="black")


#Add repeat array interval at the end
#4551649-4574359 - 1 based
tr_start=4634130-4574359
tr_end=4634130-4551649

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col="red",border="red")

#Add short stretch of repeat just after the gap

tr_start=4634130-4091213
tr_end=4634130-4087200


polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col="red",border="red")



#Add segmental duplications

start=4634130-4015195
end=4634130-4046142

polygon(x=c(start,end,end,start),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col="orange",border="orange")

start=4634130-406748
end=4634130-437745


polygon(x=c(start,end,end,start),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col="orange",border="orange")

text(y=mean(c(scaffold_north_bottom,scaffold_north_top)),x=-0.5E6,"northern",cex=0.8)

###########################################

#CHIFFCHAFF

###########################################

chr5_cc_mummer=read.delim("chr5.scaffolds.2.cc.vs.southern_ww.filt.delta.coords.out",header=F)

chr5_cc_mummer_filt=chr5_cc_mummer[which(chr5_cc_mummer[,5]>=min_aln & chr5_cc_mummer[,5]>=min_aln),]

chr5_cc_mummer_filt=chr5_cc_mummer_filt[which(chr5_cc_mummer_filt[,2]<=5E6),]

polygon(x=c(0,5E6,5E6,0),y=c(scaffold_cc_bottom,scaffold_cc_bottom,scaffold_cc_top,scaffold_cc_top),col=col_scaff2,border=col_scaff2)
text(2.5E6,scaffold_cc_bottom-0.1,"51l (+)",cex=0.7) 

#Synteny



for(i in c(1:nrow(chr5_cc_mummer_filt))){

	input=chr5_cc_mummer_filt[i,]

	target_start=as.numeric(input[1])
	target_end=as.numeric(input[2])

	query_start=as.numeric(input[3])
	query_end=as.numeric(input[4])

	polygon(x=c(query_start,query_end,target_end,target_start),c(scaffold_synteny_bottom_cc,scaffold_synteny_bottom_cc,scaffold_synteny_top_cc,scaffold_synteny_top_cc),col=synteny_col,border=synteny_col_border)
}


#Add TR position 

tr_start=244213
tr_end=287315

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_cc_bottom,scaffold_cc_bottom,scaffold_cc_top,scaffold_cc_top),col="red",border="red")

text(y=mean(c(scaffold_cc_bottom,scaffold_cc_top)),x=-0.5E6,"chiffchaff",cex=0.8)

###########################################

#FLYCATCHER

###########################################

#The flycatcher chromosome should be drawn in a reverse orientation 
#In this case, we will have to draw two separate blocks 

#The first block is from 8016461-8332205 and should be from the start
#polygon(x=c(0,(8332205-8016461),(8332205-8016461),0),y=c(scaffold_fl_bottom,scaffold_fl_bottom,scaffold_fl_top,scaffold_fl_top),col=col_scaff2,border=col_scaff2)
#text((8332205-8016461)/2+1E5,scaffold_fl_bottom-0.1,"Chr5:8.0-8.3 Mb",cex=0.7)

polygon(x=c(4.5E6,(8332205-8016461+4.5E6),(8332205-8016461+4.5E6),4.5E6),y=c(scaffold_fl_bottom,scaffold_fl_bottom,scaffold_fl_top,scaffold_fl_top),col=col_scaff2,border=col_scaff2)
text((8332205-8016461)/2+4.5E6,scaffold_fl_bottom-0.1,"Chr5:8.0-8.3 Mb (+)",cex=0.7)



#The second block is 218156-3672654. It should be drawn starting from 5E6-3672654 and ending at 5E6-218156
#polygon(x=c(5E6-3672654,5E6-218156,5E6-218156,5E6-3672654),y=c(scaffold_fl_bottom,scaffold_fl_bottom,scaffold_fl_top,scaffold_fl_top),col=col_scaff2,border=col_scaff2)
#text(5E6-mean(c(3672654,218156)),scaffold_fl_bottom-0.1,"Chr5:0.2-3.7 Mb",cex=0.7)

polygon(x=c(3.9E6-3672654,3.9E6-218156,3.9E6-218156,3.9E6-3672654),y=c(scaffold_fl_bottom,scaffold_fl_bottom,scaffold_fl_top,scaffold_fl_top),col=col_scaff2,border=col_scaff2)
text(3.9E6-mean(c(3672654,218156)),scaffold_fl_bottom-0.1,"Chr5:0.2-3.7 Mb (-)",cex=0.7)

chr5_flycatcher_mummer=read.delim("chr5.scaffolds.ww_south.vs.flycatcher.delta.filt.coords.out",header=F)

chr5_flycatcher_mummer_filt=chr5_flycatcher_mummer[which(chr5_flycatcher_mummer[,5]>=min_aln & chr5_flycatcher_mummer[,6]>=min_aln & chr5_flycatcher_mummer[,4]<=5E06),]



#Synteny

#Here we can plot the flycatcher in a reverse orientation 

#First block > 4 Mb in flycatcher...


for(i in c(1:nrow(chr5_flycatcher_mummer_filt))){

	input=chr5_flycatcher_mummer_filt[i,]

	target_start=as.numeric(input[3])
	target_end=as.numeric(input[4])

	if(as.numeric(input[1])>8E6){
		#query_start=8332205-as.numeric(input[1])+4.5E6
		#query_end=8332205-as.numeric(input[2])+4.5E6
		query_start=as.numeric(input[1])-8016461+4.5E6
		query_end=as.numeric(input[2])-8016461+4.5E6
		}else{
		#query_start=5E6-as.numeric(input[1])
		#query_end=5E6-as.numeric(input[2])
		query_start=3.9E6-as.numeric(input[1])
		query_end=3.9E6-as.numeric(input[2])
		
		}
		
	
	polygon(x=c(query_start,query_end,target_end,target_start),c(scaffold_synteny_bottom_fl,scaffold_synteny_bottom_fl,scaffold_synteny_top_fl,scaffold_synteny_top_fl),col=synteny_col,border=synteny_col_border)
}


text(y=mean(c(scaffold_fl_bottom,scaffold_fl_top)),x=-0.5E6,"flycatcher",cex=0.8)

```













# Supplementary Figure 1 

R code to generate Supplementary Figure 1 - normalized read coverage on differentiated scaffolds: 

```
high_fst_data=read.delim("fst_0.7_variants.bed",header=F)

read_counts=read.delim("ww_read_counts.southern_hifi.1kb.windows.out",header=F)
names(read_counts)=c("scaffold","start","end","1A05","UK06","0G03","0G04","0G10","0J01","1P02","3K06","7A12","96A01","96B07","1M13","1N12","1K05","1K10","1L17","1L19","1L20","1M08","1O01","1O04","1O06") 

norm_read_counts=as.data.frame(array(,dim(read_counts)))

norm_read_counts[,1:3]=read_counts[,1:3]

median_counts_samples=rep(0,22)


for(i in c(1:length(median_counts_samples))){
median_counts_samples[i]=median(read_counts[,i+3])
}


for(i in c(1:length(median_counts_samples))){
norm_read_counts[,i+3]=read_counts[,i+3]/median_counts_samples[i]
}


southern_mean_norm=rowMeans(norm_read_counts[,4:14])
northern_mean_norm=rowMeans(norm_read_counts[,15:25])

norm_read_counts=cbind(norm_read_counts,southern_mean_norm,northern_mean_norm)

scaffolds=c("Scaffold19,Scaffold94")


for(i in c(1:length(scaffolds))){

scaffold=scaffolds[i]

scaff_data=norm_read_counts[which(norm_read_counts[,1]==scaffold),]

filename=gsub(" ","",paste("coverage_plots_north_south_reseq/",scaffold,".png"))

png(filename=filename,width = 960, height = 480)

plot(rowMeans(scaff_data[,2:3]),log10(scaff_data$northern_mean_norm/scaff_data$southern_mean_norm),type="l",xlab="position",ylab="log10(norm cov north/norm cov south)",main=scaffold,ylim=c(-3,3))

high_fst_data_scaff=high_fst_data[which(high_fst_data[,1]==scaffold),]

if(nrow(high_fst_data_scaff)>0){
points(high_fst_data_scaff[,3],rep(3,nrow(high_fst_data_scaff)),col="red",pch=16)
}

dev.off()

}
```


# Supplementary Figure 2 

R code to generate Supplementary Figure 2 - linked read barcode coverage:

```
northern_chromium=read.delim("P6352_101.southern_genome.sorted.bam.molecules.1kb_coverage.bed",header=F)
southern_chromium=read.delim("P6352_102.southern_genome.sorted.bam.molecules.1kb_coverage.bed",header=F)
southern2_chromium=read.delim("P13103_302.southern_genome.sorted.bam.molecules.1kb_coverage.bed",header=F)

northern_chromium_norm=northern_chromium
southern_chromium_norm=southern_chromium
southern2_chromium_norm=southern2_chromium

northern_chromium_norm[,4]=northern_chromium[,4]/median(northern_chromium[,4])
southern_chromium_norm[,4]=southern_chromium[,4]/median(southern_chromium[,4])
southern2_chromium_norm[,4]=southern2_chromium[,4]/median(southern2_chromium[,4])



#Set a generic line width parameter
plot_lwd=2

#################################################

#Chr1

scaffold="Scaffold19"

northern_chromium_scaffold=northern_chromium_norm[which(northern_chromium_norm[,1]==scaffold),]
southern_chromium_scaffold=southern_chromium_norm[which(southern_chromium_norm[,1]==scaffold),]
southern2_chromium_scaffold=southern2_chromium_norm[which(southern2_chromium_norm[,1]==scaffold),]

#Start
plot(rowMeans(southern_chromium_scaffold[,2:3]),southern_chromium_scaffold[,4],lwd=plot_lwd,type="l",col="green",ylim=c(-0.3,1.5),xlim=c(0,1E5),xlab="position",ylab="normalized molecule coverage",main="Scaffold19 (Chr1) - start")
points(rowMeans(northern_chromium_scaffold[,2:3]),northern_chromium_scaffold[,4],lwd=plot_lwd,type="l",col="blue")
points(rowMeans(southern2_chromium_scaffold[,2:3]),southern2_chromium_scaffold[,4],lwd=plot_lwd,type="l",col="orange")

#Bed intervals of tandem repeats
#Scaffold19	0	48983	rnd-5_family-604	7553	-	14.6	2.4	3.1	(11803076)	Unknown	(0)	1598	1	224054
#Scaffold19	57733	58159	rnd-5_family-604	2344	-	13.4	0.2	1.7	(11793900)	Unknown	(808)	790	371	224058

start=1
end=48983
polygon(x=c(start,end,end,start),y=c(-0.2,-0.2,-0.1,-0.1),col="red",border="red")

start=57733+1
end=58159
polygon(x=c(start,end,end,start),y=c(-0.2,-0.2,-0.1,-0.1),col="red",border="red")

#Alignment start of northern genome based on mummer 
abline(v=58074,lty="dotted")

#Mid breakpoint region 
plot(rowMeans(southern_chromium_scaffold[,2:3]),southern_chromium_scaffold[,4],lwd=plot_lwd,type="l",col="green",ylim=c(-0.3,1.5),xlim=c(7.9E6,8.1E6),xlab="position",ylab="normalized molecule coverage",main="Scaffold19 (Chr1) - mid breakpoint region")
points(rowMeans(northern_chromium_scaffold[,2:3]),northern_chromium_scaffold[,4],lwd=plot_lwd,type="l",col="blue")
points(rowMeans(southern2_chromium_scaffold[,2:3]),southern2_chromium_scaffold[,4],lwd=plot_lwd,type="l",col="orange")

#Bed intervals of tandem repeat at this breakpoint:
#7978401 7978481
#7978792 7980340
#7980325 7980360

#First interval
start=7978401+1
end=7978481
polygon(x=c(start,end,end,start),y=c(-0.2,-0.2,-0.1,-0.1),col="red",border="red")

#Last two intervals are overlapping
start=7978792+1
end=7980360
polygon(x=c(start,end,end,start),y=c(-0.2,-0.2,-0.1,-0.1),col="red",border="red")

#Add boundary determined from genome alignment - note that the northern genome is incomplete upstream
abline(v=7968505,lty="dotted")
abline(v=7978783,lty="dotted")


#End
plot(rowMeans(southern_chromium_scaffold[,2:3]),southern_chromium_scaffold[,4],type="l",lwd=plot_lwd,col="green",ylim=c(-0.3,1.5),xlim=c(11.66E6,11.7E6),xlab="position",ylab="normalized molecule coverage",main="Scaffold19 (Chr1) - end")
points(rowMeans(northern_chromium_scaffold[,2:3]),northern_chromium_scaffold[,4],type="l",lwd=plot_lwd,col="blue")
points(rowMeans(southern2_chromium_scaffold[,2:3]),southern2_chromium_scaffold[,4],type="l",lwd=plot_lwd,col="orange")

#TR intervals (bed)
#Scaffold19	11678265	11678329	rnd-5_family-604	241	+	13.7	4.6	0.6	(173730)	Unknown	1192	1231	(367)	229602	*
#Scaffold19	11678314	11852059	rnd-5_family-604	7773	-	14.3	3.0	3.5	(0)	Unknown	(0)	1598	1	229603

#The two intervals are overlapping
start=11678265+1
end=11852059
polygon(x=c(start,end,end,start),y=c(-0.2,-0.2,-0.1,-0.1),col="red",border="red")

#Northern alignment based on mummer
abline(v=11680761,lty="dotted")






#################################################

#Chr5


scaffold="Scaffold0"

northern_chromium_scaffold=northern_chromium_norm[which(northern_chromium_norm[,1]==scaffold),]
southern_chromium_scaffold=southern_chromium_norm[which(southern_chromium_norm[,1]==scaffold),]
southern2_chromium_scaffold=southern2_chromium_norm[which(southern2_chromium_norm[,1]==scaffold),]

#Start
plot(rowMeans(southern_chromium_scaffold[,2:3]),southern_chromium_scaffold[,4],type="l",lwd=plot_lwd,col="green",ylim=c(-0.3,1.5),xlim=c(0.2E6,0.4E6),xlab="position",ylab="normalized molecule coverage",main="Scaffold0 (Chr5) - start")
points(rowMeans(northern_chromium_scaffold[,2:3]),northern_chromium_scaffold[,4],type="l",lwd=plot_lwd,col="blue")
points(rowMeans(southern2_chromium_scaffold[,2:3]),southern2_chromium_scaffold[,4],type="l",lwd=plot_lwd,col="orange")

#Position of tandem repeat array
#Scaffold0       212053  212092  rnd-6_family-6462       323     -       0.0     0.0     0.0     (66660860)      Unknown (1099)  891     853     156     *
#Scaffold0       328001  328479  rnd-6_family-6462       3731    -       0.6     0.8     0.8     (66544473)      Unknown (0)     1990    1513    253
start=212053+1
end=328479
polygon(x=c(start,end,end,start),y=c(-0.2,-0.2,-0.1,-0.1),col="red",border="red")


#Mummer boundaries
#Scaffold0	200547	203878	Scaffold68	4027155	4030464	0.950165	+ (Satsuma)
#200402	210250	418869	428707	9849	9839	92.06	1	1	Scaffold0	Scaffold68
#211396	212092	429030	429697	697	668	92.68	1	1	Scaffold0	Scaffold68
#327097	467392	430047	570373	140296	140327	98.19	1	1	Scaffold0	Scaffold68
abline(v=200402,lty="dotted")
abline(v=212092,lty="dotted")
abline(v=327097,lty="dotted")

#Mid part
plot(rowMeans(southern_chromium_scaffold[,2:3]),southern_chromium_scaffold[,4],type="l",lwd=plot_lwd,col="green",ylim=c(0,1.5),xlim=c(3.8E6,4E6),xlab="position",ylab="normalized molecule coverage",main="Scaffold0 (Chr5) - mid part")
points(rowMeans(northern_chromium_scaffold[,2:3]),northern_chromium_scaffold[,4],type="l",lwd=plot_lwd,col="blue")
points(rowMeans(southern2_chromium_scaffold[,2:3]),southern2_chromium_scaffold[,4],type="l",lwd=plot_lwd,col="orange")
#Mummer boundaries
#3905830	3915477	4000447	4010209	9648	9763	95.09	1	1	Scaffold0	Scaffold68
#3917392	3926506	4541627	4532566	9115	9062	95.85	1	-1	Scaffold0	Scaffold68
abline(v=3905830,lty="dotted")
abline(v=3917392,lty="dotted")

#End part
plot(rowMeans(southern_chromium_scaffold[,2:3]),southern_chromium_scaffold[,4],type="l",lwd=plot_lwd,col="green",ylim=c(0,1.5),xlim=c(4.3E6,4.4E6),xlab="position",ylab="normalized molecule coverage",main="Scaffold0 (Chr5) - end part")
points(rowMeans(northern_chromium_scaffold[,2:3]),northern_chromium_scaffold[,4],type="l",lwd=plot_lwd,col="blue")
points(rowMeans(southern2_chromium_scaffold[,2:3]),southern2_chromium_scaffold[,4],type="l",lwd=plot_lwd,col="orange")
#Mummer boundary
#4358697	4361809	4101323	4098206	3113	3118	95.88	1	-1	Scaffold0	Scaffold68
abline(v=4361809,lty="dotted")
```


