# Figure 2 

R code to generate Figure 2:

```
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

col_scaff2="seagreen2"
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
scaffold_synteny_bottom_ww_north=-0.9+0.02

scaffold_synteny_top_cc=-1.2
scaffold_synteny_bottom_cc=-1.5+0.02


scaffold_synteny_top_fl=-1.8
scaffold_synteny_bottom_fl=-2.1+0.02


synteny_col="lightgoldenrod"
synteny_col_border="NA"
synteny_col_border=synteny_col

min_aln=2000

#Repeat colors
rep_color="blue"
dupl_interval_color="#999999"




in2mm=25.4

pdf("Figure2.pdf",width=180/in2mm,height=185/in2mm)


layout(matrix(1:3, ncol = 1), widths = 1, heights = c(0.5,0.5,0.5,0.5,0.6), respect = FALSE)
par(mar = c(0.3, 4.1, 2, 2.1))


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

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col=rep_color,border=rep_color)

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
polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col=rep_color,border=rep_color)
#End:11678266-11852059 ~ 174 kb - 1-based
tr_start=11852059-11852059+cumulative_pos
tr_end=11852059-11678265+cumulative_pos
polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col=rep_color,border=rep_color)
#Middle: 7978793-7980340 - 1-based (1548 bp)
tr_start=11852059-7980340+cumulative_pos
tr_end=11852059-7978793+cumulative_pos
polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col=rep_color,border=rep_color)

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
polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col=rep_color,border=rep_color)

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

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col=rep_color,border=rep_color)

#TR in the middle - here we can include the gap (1-based)
#Scaffold156     3675482 3685135
#Scaffold156     3742954 3743449 

tr_start=cumulative_pos+3675482
tr_end=cumulative_pos+3743449

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col=rep_color,border=rep_color)

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

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_cc_bottom,scaffold_cc_bottom,scaffold_cc_top,scaffold_cc_top),col=rep_color,border=rep_color)


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

#Plot scaffold segment 
polygon(x=c(cumulative_pos,cumulative_pos+14E6,cumulative_pos+14E6,cumulative_pos),y=c(scaffold_cc_bottom,scaffold_cc_bottom,scaffold_cc_top,scaffold_cc_top),col=col_scaff2,border=col_scaff2)
text(mean(c(cumulative_pos,cumulative_pos+14E6)),scaffold_cc_bottom-0.1,"7l (-)",cex=0.7)

#Add TR at the end of the div region: 6512119-6952665 (1-based) ~440 kb

tr_start=cumulative_pos+19027237-6512120
tr_end=cumulative_pos+19027237-6952665

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_cc_bottom,scaffold_cc_bottom,scaffold_cc_top,scaffold_cc_top),col=rep_color,border=rep_color)

#Add TR in the central part of the div region: 14900329-15170453 (1-based) ~ 270 kb

tr_start=cumulative_pos+19027237-14900328
tr_end=cumulative_pos+19027237-15170453

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_cc_bottom,scaffold_cc_bottom,scaffold_cc_top,scaffold_cc_top),col=rep_color,border=rep_color)


#Add TR at the start of the div region. Here there are two close intervals: 
#1) 18874730 - 18983007 (1 based) ~ 108 kb

tr_start=cumulative_pos+19027237-18874730
tr_end=cumulative_pos+19027237-18983007

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_cc_bottom,scaffold_cc_bottom,scaffold_cc_top,scaffold_cc_top),col=rep_color,border=rep_color)

#2 18989765 - 19027234 (1 based) - 37.5 kb 

tr_start=cumulative_pos+19027237-18989765
tr_end=cumulative_pos+19027237-19027234

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_cc_bottom,scaffold_cc_bottom,scaffold_cc_top,scaffold_cc_top),col=rep_color,border=rep_color)



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
text(plot_end/2,scaffold_fl_bottom-0.1,"Chr1:62.5-78.6 Mb (+)",cex=0.7)

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

text(y=0.5,x=-1.9E6,"FST",srt=90,cex=0.8)


#Plot scaffold

start_scaffold=0
end_scaffold=interval_size
polygon(x=c(start_scaffold,end_scaffold,end_scaffold,start_scaffold),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col=col_scaff1,border=col_scaff1)
text(1E6,scaffold_bottom-0.1,"38 (+)",cex=0.7)

#Plot TR - 35897021-35991244 (1-based) ~ 94224

tr_start=35897021-(35991252-interval_size)
tr_end=35991244-(35991252-interval_size)

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col=rep_color,border=rep_color)

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
polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col=rep_color,border=rep_color)

#2) Scaffold61:69204988-69283336 (1-based) ~ 78349 bp
tr_start=69302477-69204988+cumulative_pos
tr_end=69302477-69283336+cumulative_pos

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col=rep_color,border=rep_color)


#TR at the other end of the divergent region: 55894452-56078109 (1-based) ~ 183,658 bp
tr_start=69302477-55894451+cumulative_pos
tr_end=69302477-56078109+cumulative_pos

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col=rep_color,border=rep_color)


#Add line to show the divergent region - we could define it based on the synteny
#56069986-69205658
segments(69302477-69205658+cumulative_pos,1.1,69302477-56069986+cumulative_pos,1.1)
text(mean(c(69302477-69205658,69302477-56069986))+cumulative_pos,1.18,"13.1 Mb",cex=0.7)


###########################################

#NORTHERN WW

###########################################

scaffold_synteny_top=-0.4
scaffold_synteny_bottom=-0.65

##############

#Upstream scaffold - Scaffold139 (reverse) - 32,231,478 bp

start_scaffold=0
end_scaffold=interval_size

polygon(x=c(start_scaffold,end_scaffold,end_scaffold,start_scaffold),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col=col_scaff1,border=col_scaff1)
text(1E6,scaffold_north_bottom-0.1,"139 (-)",cex=0.7)


#Add TR - 1-3135 (1-based)

tr_start=interval_size-3135
tr_end=interval_size

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col=rep_color,border=rep_color)

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

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col=rep_color,border=rep_color)



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

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col=rep_color,border=rep_color)


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

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_cc_bottom,scaffold_cc_bottom,scaffold_cc_top,scaffold_cc_top),col=rep_color,border=rep_color)


#ptg000040l:29049553-29328563 (1-based) - 279,011 bp

tr_start=29049552-(29328574-15.5E6)
tr_end=29328563-(29328574-15.5E6)

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_cc_bottom,scaffold_cc_bottom,scaffold_cc_top,scaffold_cc_top),col=rep_color,border=rep_color)

#Final TR at the very end: 29269249-29328563 (1-based) - 59,315 bp

tr_start=29269249-(29328574-15.5E6)
tr_end=29328563-(29328574-15.5E6)

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_cc_bottom,scaffold_cc_bottom,scaffold_cc_top,scaffold_cc_top),col=rep_color,border=rep_color)



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

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_cc_bottom,scaffold_cc_bottom,scaffold_cc_top,scaffold_cc_top),col=rep_color,border=rep_color)


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


polygon(x=c(0,plot_end,plot_end,0),y=c(scaffold_fl_bottom,scaffold_fl_bottom,scaffold_fl_top,scaffold_fl_top),col=col_scaff2,border=col_scaff2)
text(plot_end/2,scaffold_fl_bottom-0.1,"Chr3:42.4-60.0 Mb (+)",cex=0.7)

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
#axis(2,pos=-2E5,at=c(0:10)/10,cex.axis=0.7)
axis(2,pos=-3.5E5,at=c(0:10)/10,cex.axis=0.7)

text(y=0.5,x=-0.8E6,"FST",srt=90,cex=0.8)

#Plot fst for 10 kb windows
points(rowMeans(fst_10kb_maf_scaff[,2:3]),fst_10kb_maf_scaff[,5],type="l",col=fst_10kb_plot_col,lwd=fst_10kb_plot_lwd)


#Plot scaffold segment 

polygon(x=c(0,plot_end,plot_end,0),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col=col_scaff2,border=col_scaff2)
text(mean(c(plot_start,plot_end)),scaffold_bottom-0.1,"0 (+)",cex=0.7)


#TR at the start 212054-328479 (1-based) - 116 kb

tr_start=212054
tr_end=328479

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_bottom,scaffold_bottom,scaffold_top,scaffold_top),col=rep_color,border=rep_color)


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

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col=rep_color,border=rep_color)

#Add short stretch of repeat just after the gap

tr_start=4634130-4091213
tr_end=4634130-4087200


polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col=rep_color,border=rep_color)



#Add segmental duplications

start=4634130-4015195
end=4634130-4046142

polygon(x=c(start,end,end,start),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col=dupl_interval_color,border=dupl_interval_color)

start=4634130-406748
end=4634130-437745


polygon(x=c(start,end,end,start),y=c(scaffold_north_bottom,scaffold_north_bottom,scaffold_north_top,scaffold_north_top),col=dupl_interval_color,border=dupl_interval_color)

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

polygon(x=c(tr_start,tr_end,tr_end,tr_start),y=c(scaffold_cc_bottom,scaffold_cc_bottom,scaffold_cc_top,scaffold_cc_top),col=rep_color,border=rep_color)

text(y=mean(c(scaffold_cc_bottom,scaffold_cc_top)),x=-0.5E6,"chiffchaff",cex=0.8)

###########################################

#FLYCATCHER

###########################################

#The flycatcher chromosome should be drawn in a reverse orientation 
#In this case, we will have to draw two separate blocks 

#The first block is from 8016461-8332205 and should be from the start
polygon(x=c(4.5E6,(8332205-8016461+4.5E6),(8332205-8016461+4.5E6),4.5E6),y=c(scaffold_fl_bottom,scaffold_fl_bottom,scaffold_fl_top,scaffold_fl_top),col=col_scaff2,border=col_scaff2)
text((8332205-8016461)/2+4.5E6,scaffold_fl_bottom-0.1,"Chr5:8.0-8.3 Mb (+)",cex=0.7)



#The second block is 218156-3672654. It should be drawn starting from 5E6-3672654 and ending at 5E6-218156
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
		query_start=as.numeric(input[1])-8016461+4.5E6
		query_end=as.numeric(input[2])-8016461+4.5E6
		}else{
		query_start=3.9E6-as.numeric(input[1])
		query_end=3.9E6-as.numeric(input[2])
		
		}
		
	
	polygon(x=c(query_start,query_end,target_end,target_start),c(scaffold_synteny_bottom_fl,scaffold_synteny_bottom_fl,scaffold_synteny_top_fl,scaffold_synteny_top_fl),col=synteny_col,border=synteny_col_border)
}


text(y=mean(c(scaffold_fl_bottom,scaffold_fl_top)),x=-0.5E6,"flycatcher",cex=0.8)

dev.off()
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

#Set plot colors

northern_col="blue"
southern_col="green"
southern2_col="#AA4499"
repeat_color="dark blue"


plotdir="figures/"


#################################################

#Chr1

scaffold="Scaffold19"

northern_chromium_scaffold=northern_chromium_norm[which(northern_chromium_norm[,1]==scaffold),]
southern_chromium_scaffold=southern_chromium_norm[which(southern_chromium_norm[,1]==scaffold),]
southern2_chromium_scaffold=southern2_chromium_norm[which(southern2_chromium_norm[,1]==scaffold),]

#Start
pdf(gsub(" ","",paste(plotdir,"FigureS2.1.pdf")),width=10)
plot(rowMeans(southern_chromium_scaffold[,2:3]),southern_chromium_scaffold[,4],lwd=plot_lwd,type="l",col=southern_col,ylim=c(-0.3,1.5),xlim=c(0,1E5),xlab="position",ylab="normalized molecule coverage",main="Scaffold19 (Chr1) - start")
points(rowMeans(northern_chromium_scaffold[,2:3]),northern_chromium_scaffold[,4],lwd=plot_lwd,type="l",col=northern_col)
points(rowMeans(southern2_chromium_scaffold[,2:3]),southern2_chromium_scaffold[,4],lwd=plot_lwd,type="l",col=southern2_col)

#Bed intervals of tandem repeats
#Scaffold19	0	48983	rnd-5_family-604	7553	-	14.6	2.4	3.1	(11803076)	Unknown	(0)	1598	1	224054
#Scaffold19	57733	58159	rnd-5_family-604	2344	-	13.4	0.2	1.7	(11793900)	Unknown	(808)	790	371	224058

start=1
end=48983
polygon(x=c(start,end,end,start),y=c(-0.2,-0.2,-0.1,-0.1),col=repeat_color,border=repeat_color)

start=57733+1
end=58159
polygon(x=c(start,end,end,start),y=c(-0.2,-0.2,-0.1,-0.1),col=repeat_color,border=repeat_color)

#Alignment start of northern genome based on mummer 
abline(v=58074,lty="dotted",lwd=2)
dev.off()

##########
#Mid breakpoint region

pdf(gsub(" ","",paste(plotdir,"FigureS2.2.pdf")),width=10) 
plot(rowMeans(southern_chromium_scaffold[,2:3]),southern_chromium_scaffold[,4],lwd=plot_lwd,type="l",col=southern_col,ylim=c(-0.3,1.5),xlim=c(7.9E6,8.1E6),xlab="position",ylab="normalized molecule coverage",main="Scaffold19 (Chr1) - mid breakpoint region")
points(rowMeans(northern_chromium_scaffold[,2:3]),northern_chromium_scaffold[,4],lwd=plot_lwd,type="l",col=northern_col)
points(rowMeans(southern2_chromium_scaffold[,2:3]),southern2_chromium_scaffold[,4],lwd=plot_lwd,type="l",col=southern2_col)

#Bed intervals of tandem repeat at this breakpoint:
#7978401 7978481
#7978792 7980340
#7980325 7980360

#First interval
start=7978401+1
end=7978481
polygon(x=c(start,end,end,start),y=c(-0.2,-0.2,-0.1,-0.1),col=repeat_color,border=repeat_color)

#Last two intervals are overlapping
start=7978792+1
end=7980360
polygon(x=c(start,end,end,start),y=c(-0.2,-0.2,-0.1,-0.1),col=repeat_color,border=repeat_color)

#Add boundary determined from genome alignment - note that the northern genome is incomplete upstream
abline(v=7968505,lty="dotted",lwd=2)
abline(v=7978783,lty="dotted",lwd=2)
dev.off()


##########
#End
pdf(gsub(" ","",paste(plotdir,"FigureS2.3.pdf")),width=10)
plot(rowMeans(southern_chromium_scaffold[,2:3]),southern_chromium_scaffold[,4],type="l",lwd=plot_lwd,col=southern_col,ylim=c(-0.3,1.5),xlim=c(11.66E6,11.7E6),xlab="position",ylab="normalized molecule coverage",main="Scaffold19 (Chr1) - end")
points(rowMeans(northern_chromium_scaffold[,2:3]),northern_chromium_scaffold[,4],type="l",lwd=plot_lwd,col=northern_col)
points(rowMeans(southern2_chromium_scaffold[,2:3]),southern2_chromium_scaffold[,4],type="l",lwd=plot_lwd,col=southern2_col)

#TR intervals (bed)
#Scaffold19	11678265	11678329	rnd-5_family-604	241	+	13.7	4.6	0.6	(173730)	Unknown	1192	1231	(367)	229602	*
#Scaffold19	11678314	11852059	rnd-5_family-604	7773	-	14.3	3.0	3.5	(0)	Unknown	(0)	1598	1	229603

#The two intervals are overlapping
start=11678265+1
end=11852059
polygon(x=c(start,end,end,start),y=c(-0.2,-0.2,-0.1,-0.1),col=repeat_color,border=repeat_color)

#Northern alignment based on mummer
abline(v=11680761,lty="dotted",lwd=2)
dev.off()





#################################################

#Chr5


scaffold="Scaffold0"

northern_chromium_scaffold=northern_chromium_norm[which(northern_chromium_norm[,1]==scaffold),]
southern_chromium_scaffold=southern_chromium_norm[which(southern_chromium_norm[,1]==scaffold),]
southern2_chromium_scaffold=southern2_chromium_norm[which(southern2_chromium_norm[,1]==scaffold),]

#Start
pdf(gsub(" ","",paste(plotdir,"FigureS2.4.pdf")),width=10)
plot(rowMeans(southern_chromium_scaffold[,2:3]),southern_chromium_scaffold[,4],type="l",lwd=plot_lwd,col=southern_col,ylim=c(-0.3,1.5),xlim=c(0.2E6,0.4E6),xlab="position",ylab="normalized molecule coverage",main="Scaffold0 (Chr5) - start")
points(rowMeans(northern_chromium_scaffold[,2:3]),northern_chromium_scaffold[,4],type="l",lwd=plot_lwd,col=northern_col)
points(rowMeans(southern2_chromium_scaffold[,2:3]),southern2_chromium_scaffold[,4],type="l",lwd=plot_lwd,col=southern2_col)

#Position of tandem repeat array
#Scaffold0       212053  212092  rnd-6_family-6462       323     -       0.0     0.0     0.0     (66660860)      Unknown (1099)  891     853     156     *
#Scaffold0       328001  328479  rnd-6_family-6462       3731    -       0.6     0.8     0.8     (66544473)      Unknown (0)     1990    1513    253
start=212053+1
end=328479
polygon(x=c(start,end,end,start),y=c(-0.2,-0.2,-0.1,-0.1),col=repeat_color,border=repeat_color)


#Mummer boundaries
#Scaffold0	200547	203878	Scaffold68	4027155	4030464	0.950165	+ (Satsuma)
#200402	210250	418869	428707	9849	9839	92.06	1	1	Scaffold0	Scaffold68
#211396	212092	429030	429697	697	668	92.68	1	1	Scaffold0	Scaffold68
#327097	467392	430047	570373	140296	140327	98.19	1	1	Scaffold0	Scaffold68
abline(v=200402,lty="dotted",lwd=2)
abline(v=212092,lty="dotted",lwd=2)
abline(v=327097,lty="dotted",lwd=2)
dev.off()



#Mid part
pdf(gsub(" ","",paste(plotdir,"FigureS2.5.pdf")),width=10)
plot(rowMeans(southern_chromium_scaffold[,2:3]),southern_chromium_scaffold[,4],type="l",lwd=plot_lwd,col=southern_col,ylim=c(0,1.5),xlim=c(3.8E6,4E6),xlab="position",ylab="normalized molecule coverage",main="Scaffold0 (Chr5) - mid part")
points(rowMeans(northern_chromium_scaffold[,2:3]),northern_chromium_scaffold[,4],type="l",lwd=plot_lwd,col=northern_col)
points(rowMeans(southern2_chromium_scaffold[,2:3]),southern2_chromium_scaffold[,4],type="l",lwd=plot_lwd,col=southern2_col)
#Mummer boundaries
#3905830	3915477	4000447	4010209	9648	9763	95.09	1	1	Scaffold0	Scaffold68
#3917392	3926506	4541627	4532566	9115	9062	95.85	1	-1	Scaffold0	Scaffold68
abline(v=3905830,lty="dotted",lwd=2)
abline(v=3917392,lty="dotted",lwd=2)
dev.off()

#End part
pdf(gsub(" ","",paste(plotdir,"FigureS2.6.pdf")),width=10)
plot(rowMeans(southern_chromium_scaffold[,2:3]),southern_chromium_scaffold[,4],type="l",lwd=plot_lwd,col=southern_col,ylim=c(0,1.5),xlim=c(4.3E6,4.4E6),xlab="position",ylab="normalized molecule coverage",main="Scaffold0 (Chr5) - end part")
points(rowMeans(northern_chromium_scaffold[,2:3]),northern_chromium_scaffold[,4],type="l",lwd=plot_lwd,col=northern_col)
points(rowMeans(southern2_chromium_scaffold[,2:3]),southern2_chromium_scaffold[,4],type="l",lwd=plot_lwd,col=southern2_col)
#Mummer boundary
#4358697	4361809	4101323	4098206	3113	3118	95.88	1	-1	Scaffold0	Scaffold68
abline(v=4361809,lty="dotted",lwd=2)
dev.off()
```


# Supplementary Figure 7

R code to generate Supplementary Figure 7:

```

library(ggplot2)
library(gridExtra)


#Read the datsets

north_tajd=read.delim("freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.tajd_10kb.out",header=T)
south_tajd=read.delim("freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.tajd_10kb.out",header=T)

south_pi=read.delim("freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.pi_10kb.out",header=T)
north_pi=read.delim("freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.pi_10kb.out",header=T)

north_ld=read.delim("ww_southern_hifi_bionano.filt.northern.ld.out",header=F)
south_ld_subset=read.delim("ww_southern_hifi_bionano.filt.southern.7_samples.ld.out",header=F)

#Fay and Wu's H
south_popgenome=read.delim("southern_with_anc_popgenome.2.out",header=F)
north_popgenome=read.delim("northern_with_anc_popgenome.2.out",header=F)


###############################
#Define boundaries for chromosomes 3 and 5
chrom3_start=56E6
chrom5_start=3.27E5
chrom5_end=4.36E6




#Tajima's D

#North
north_tajd_chr1=north_tajd[which(north_tajd$CHROM=="Scaffold19"),]
north_tajd_rest=north_tajd[-which(north_tajd$CHROM=="Scaffold19"),]
north_tajd_chr3=north_tajd_rest[which(north_tajd_rest$CHROM=="Scaffold61" & north_tajd_rest$BIN_START>chrom3_start),]
north_tajd_rest=north_tajd_rest[-which(north_tajd_rest$CHROM=="Scaffold61" & north_tajd_rest$BIN_START>chrom3_start),]

north_tajd_chr5=north_tajd_rest[which(north_tajd_rest$CHROM=="Scaffold0" & north_tajd_rest$BIN_START<chrom5_end & north_tajd_rest$BIN_START>chrom5_start),]
north_tajd_rest=north_tajd_rest[-which(north_tajd_rest$CHROM=="Scaffold0" & north_tajd_rest$BIN_START<chrom5_end & north_tajd_rest$BIN_START>chrom5_start),]

north_tajd_chr1=cbind(rep("chr1",nrow(north_tajd_chr1)),north_tajd_chr1)
names(north_tajd_chr1)[1]="region"
north_tajd_chr3=cbind(rep("chr3",nrow(north_tajd_chr3)),north_tajd_chr3)
names(north_tajd_chr3)[1]="region"
north_tajd_chr5=cbind(rep("chr5",nrow(north_tajd_chr5)),north_tajd_chr5)
names(north_tajd_chr5)[1]="region"
north_tajd_rest=cbind(rep("rest of genome",nrow(north_tajd_rest)),north_tajd_rest)
names(north_tajd_rest)[1]="region"
north_tajd_comb=rbind(north_tajd_chr1,north_tajd_chr3,north_tajd_chr5,north_tajd_rest)


#South
south_tajd_chr1=south_tajd[which(south_tajd$CHROM=="Scaffold19"),]
south_tajd_rest=south_tajd[-which(south_tajd$CHROM=="Scaffold19"),]
south_tajd_chr3=south_tajd_rest[which(south_tajd_rest$CHROM=="Scaffold61" & south_tajd_rest$BIN_START>chrom3_start),]
south_tajd_rest=south_tajd_rest[-which(south_tajd_rest$CHROM=="Scaffold61" & south_tajd_rest$BIN_START>chrom3_start),]

south_tajd_chr5=south_tajd_rest[which(south_tajd_rest$CHROM=="Scaffold0" & south_tajd_rest$BIN_START<chrom5_end & south_tajd_rest$BIN_START>chrom5_start),]
south_tajd_rest=south_tajd_rest[-which(south_tajd_rest$CHROM=="Scaffold0" & south_tajd_rest$BIN_START<chrom5_end & south_tajd_rest$BIN_START>chrom5_start),]

south_tajd_chr1=cbind(rep("chr1",nrow(south_tajd_chr1)),south_tajd_chr1)
names(south_tajd_chr1)[1]="region"
south_tajd_chr3=cbind(rep("chr3",nrow(south_tajd_chr3)),south_tajd_chr3)
names(south_tajd_chr3)[1]="region"
south_tajd_chr5=cbind(rep("chr5",nrow(south_tajd_chr5)),south_tajd_chr5)
names(south_tajd_chr5)[1]="region"
south_tajd_rest=cbind(rep("rest of genome",nrow(south_tajd_rest)),south_tajd_rest)
names(south_tajd_rest)[1]="region"
south_tajd_comb=rbind(south_tajd_chr1,south_tajd_chr3,south_tajd_chr5,south_tajd_rest)



#Nucleotide diversity

#North
north_pi_chr1=north_pi[which(north_pi$CHROM=="Scaffold19"),]
north_pi_rest=north_pi[-which(north_pi$CHROM=="Scaffold19"),]
north_pi_chr3=north_pi_rest[which(north_pi_rest$CHROM=="Scaffold61" & north_pi_rest$BIN_START>chrom3_start),]
north_pi_rest=north_pi_rest[-which(north_pi_rest$CHROM=="Scaffold61" & north_pi_rest$BIN_START>chrom3_start),]

north_pi_chr5=north_pi_rest[which(north_pi_rest$CHROM=="Scaffold0" & north_pi_rest$BIN_START<chrom5_end & north_pi_rest$BIN_START>chrom5_start),]
north_pi_rest=north_pi_rest[-which(north_pi_rest$CHROM=="Scaffold0" & north_pi_rest$BIN_START<chrom5_end & north_pi_rest$BIN_START>chrom5_start),]

north_pi_chr1=cbind(rep("chr1",nrow(north_pi_chr1)),north_pi_chr1)
names(north_pi_chr1)[1]="region"
north_pi_chr3=cbind(rep("chr3",nrow(north_pi_chr3)),north_pi_chr3)
names(north_pi_chr3)[1]="region"
north_pi_chr5=cbind(rep("chr5",nrow(north_pi_chr5)),north_pi_chr5)
names(north_pi_chr5)[1]="region"
north_pi_rest=cbind(rep("rest of genome",nrow(north_pi_rest)),north_pi_rest)
names(north_pi_rest)[1]="region"
north_pi_comb=rbind(north_pi_chr1,north_pi_chr3,north_pi_chr5,north_pi_rest)


#South
south_pi_chr1=south_pi[which(south_pi$CHROM=="Scaffold19"),]
south_pi_rest=south_pi[-which(south_pi$CHROM=="Scaffold19"),]
south_pi_chr3=south_pi_rest[which(south_pi_rest$CHROM=="Scaffold61" & south_pi_rest$BIN_START>chrom3_start),]
south_pi_rest=south_pi_rest[-which(south_pi_rest$CHROM=="Scaffold61" & south_pi_rest$BIN_START>chrom3_start),]

south_pi_chr5=south_pi_rest[which(south_pi_rest$CHROM=="Scaffold0" & south_pi_rest$BIN_START<chrom5_end & south_pi_rest$BIN_START>chrom5_start),]
south_pi_rest=south_pi_rest[-which(south_pi_rest$CHROM=="Scaffold0" & south_pi_rest$BIN_START<chrom5_end & south_pi_rest$BIN_START>chrom5_start),]

south_pi_chr1=cbind(rep("chr1",nrow(south_pi_chr1)),south_pi_chr1)
names(south_pi_chr1)[1]="region"
south_pi_chr3=cbind(rep("chr3",nrow(south_pi_chr3)),south_pi_chr3)
names(south_pi_chr3)[1]="region"
south_pi_chr5=cbind(rep("chr5",nrow(south_pi_chr5)),south_pi_chr5)
names(south_pi_chr5)[1]="region"
south_pi_rest=cbind(rep("rest of genome",nrow(south_pi_rest)),south_pi_rest)
names(south_pi_rest)[1]="region"
south_pi_comb=rbind(south_pi_chr1,south_pi_chr3,south_pi_chr5,south_pi_rest)




#Linkage disequilibrium

#North
north_ld_chr1=north_ld[which(north_ld[,1]=="Scaffold19"),]
north_ld_rest=north_ld[-which(north_ld[,1]=="Scaffold19"),]
north_ld_chr3=north_ld_rest[which(north_ld_rest[,1]=="Scaffold61" & north_ld_rest[,2]>chrom3_start),]
north_ld_rest=north_ld_rest[-which(north_ld_rest[,1]=="Scaffold61" & north_ld_rest[,2]>chrom3_start),]

north_ld_chr5=north_ld_rest[which(north_ld_rest[,1]=="Scaffold0" & north_ld_rest[,2]<chrom5_end & north_ld_rest[,2]>chrom5_start),]
north_ld_rest=north_ld_rest[-which(north_ld_rest[,1]=="Scaffold0" & north_ld_rest[,2]<chrom5_end & north_ld_rest[,2]>chrom5_start),]

north_ld_chr1=cbind(rep("chr1",nrow(north_ld_chr1)),north_ld_chr1)
names(north_ld_chr1)[1]="region"
north_ld_chr3=cbind(rep("chr3",nrow(north_ld_chr3)),north_ld_chr3)
names(north_ld_chr3)[1]="region"
north_ld_chr5=cbind(rep("chr5",nrow(north_ld_chr5)),north_ld_chr5)
names(north_ld_chr5)[1]="region"
north_ld_rest=cbind(rep("rest of genome",nrow(north_ld_rest)),north_ld_rest)
names(north_ld_rest)[1]="region"
north_ld_comb=rbind(north_ld_chr1,north_ld_chr3,north_ld_chr5,north_ld_rest)

names(north_ld_comb)=c("region","scaffold","start","end","r2","Dprime")


#South subset
south_ld_subset_chr1=south_ld_subset[which(south_ld_subset[,1]=="Scaffold19"),]
south_ld_subset_rest=south_ld_subset[-which(south_ld_subset[,1]=="Scaffold19"),]
south_ld_subset_chr3=south_ld_subset_rest[which(south_ld_subset_rest[,1]=="Scaffold61" & south_ld_subset_rest[,2]>chrom3_start),]
south_ld_subset_rest=south_ld_subset_rest[-which(south_ld_subset_rest[,1]=="Scaffold61" & south_ld_subset_rest[,2]>chrom3_start),]

south_ld_subset_chr5=south_ld_subset_rest[which(south_ld_subset_rest[,1]=="Scaffold0" & south_ld_subset_rest[,2]<chrom5_end & south_ld_subset_rest[,2]>chrom5_start),]
south_ld_subset_rest=south_ld_subset_rest[-which(south_ld_subset_rest[,1]=="Scaffold0" & south_ld_subset_rest[,2]<chrom5_end & south_ld_subset_rest[,2]>chrom5_start),]

south_ld_subset_chr1=cbind(rep("chr1",nrow(south_ld_subset_chr1)),south_ld_subset_chr1)
names(south_ld_subset_chr1)[1]="region"
south_ld_subset_chr3=cbind(rep("chr3",nrow(south_ld_subset_chr3)),south_ld_subset_chr3)
names(south_ld_subset_chr3)[1]="region"
south_ld_subset_chr5=cbind(rep("chr5",nrow(south_ld_subset_chr5)),south_ld_subset_chr5)
names(south_ld_subset_chr5)[1]="region"
south_ld_subset_rest=cbind(rep("rest of genome",nrow(south_ld_subset_rest)),south_ld_subset_rest)
names(south_ld_subset_rest)[1]="region"
south_ld_subset_comb=rbind(south_ld_subset_chr1,south_ld_subset_chr3,south_ld_subset_chr5,south_ld_subset_rest)

names(south_ld_subset_comb)=c("region","scaffold","start","end","r2","Dprime")


#Fay and Wu's H

#South
south_popgenome_chr1=south_popgenome[which(south_popgenome[,1]=="Scaffold19"),]
south_popgenome_rest=south_popgenome[-which(south_popgenome[,1]=="Scaffold19"),]
south_popgenome_chr3=south_popgenome_rest[which(south_popgenome_rest[,1]=="Scaffold61" & south_popgenome_rest[,2]>chrom3_start),]
south_popgenome_rest=south_popgenome_rest[-which(south_popgenome_rest[,1]=="Scaffold61" & south_popgenome_rest[,2]>chrom3_start),]

south_popgenome_chr5=south_popgenome_rest[which(south_popgenome_rest[,1]=="Scaffold0" & south_popgenome_rest[,2]<chrom5_end & south_popgenome_rest[,2]>chrom5_start),]
south_popgenome_rest=south_popgenome_rest[-which(south_popgenome_rest[,1]=="Scaffold0" & south_popgenome_rest[,2]<chrom5_end & south_popgenome_rest[,2]>chrom5_start),]

south_popgenome_chr1=cbind(rep("chr1",nrow(south_popgenome_chr1)),south_popgenome_chr1)
names(south_popgenome_chr1)[1]="region"
south_popgenome_chr3=cbind(rep("chr3",nrow(south_popgenome_chr3)),south_popgenome_chr3)
names(south_popgenome_chr3)[1]="region"
south_popgenome_chr5=cbind(rep("chr5",nrow(south_popgenome_chr5)),south_popgenome_chr5)
names(south_popgenome_chr5)[1]="region"
south_popgenome_rest=cbind(rep("rest of genome",nrow(south_popgenome_rest)),south_popgenome_rest)
names(south_popgenome_rest)[1]="region"
south_popgenome_comb=rbind(south_popgenome_chr1,south_popgenome_chr3,south_popgenome_chr5,south_popgenome_rest)

names(south_popgenome_comb)=c("region","Scaffold","start","end","tajd","H")

#North
north_popgenome_chr1=north_popgenome[which(north_popgenome[,1]=="Scaffold19"),]
north_popgenome_rest=north_popgenome[-which(north_popgenome[,1]=="Scaffold19"),]
north_popgenome_chr3=north_popgenome_rest[which(north_popgenome_rest[,1]=="Scaffold61" & north_popgenome_rest[,2]>chrom3_start),]
north_popgenome_rest=north_popgenome_rest[-which(north_popgenome_rest[,1]=="Scaffold61" & north_popgenome_rest[,2]>chrom3_start),]

north_popgenome_chr5=north_popgenome_rest[which(north_popgenome_rest[,1]=="Scaffold0" & north_popgenome_rest[,2]<chrom5_end & north_popgenome_rest[,2]>chrom5_start),]
north_popgenome_rest=north_popgenome_rest[-which(north_popgenome_rest[,1]=="Scaffold0" & north_popgenome_rest[,2]<chrom5_end & north_popgenome_rest[,2]>chrom5_start),]

north_popgenome_chr1=cbind(rep("chr1",nrow(north_popgenome_chr1)),north_popgenome_chr1)
names(north_popgenome_chr1)[1]="region"
north_popgenome_chr3=cbind(rep("chr3",nrow(north_popgenome_chr3)),north_popgenome_chr3)
names(north_popgenome_chr3)[1]="region"
north_popgenome_chr5=cbind(rep("chr5",nrow(north_popgenome_chr5)),north_popgenome_chr5)
names(north_popgenome_chr5)[1]="region"
north_popgenome_rest=cbind(rep("rest of genome",nrow(north_popgenome_rest)),north_popgenome_rest)
names(north_popgenome_rest)[1]="region"
north_popgenome_comb=rbind(north_popgenome_chr1,north_popgenome_chr3,north_popgenome_chr5,north_popgenome_rest)

names(north_popgenome_comb)=c("region","Scaffold","start","end","tajd","H")


########################################

#Create the different plots



tajd_north=ggplot(north_tajd_comb, aes(x=TajimaD, fill=region)) + geom_density(alpha=.3) + theme_classic() + xlim(-3,3) + ggtitle("Tajima's D northern") + theme(plot.title = element_text(hjust = 0.5)) + labs(x = "D")
tajd_south=ggplot(south_tajd_comb, aes(x=TajimaD, fill=region)) + geom_density(alpha=.3) + theme_classic() + xlim(-3,3) + ggtitle("Tajima's D southern") + theme(plot.title = element_text(hjust = 0.5)) + labs(x = "D")

pi_north=ggplot(north_pi_comb, aes(x=PI, fill=region)) + geom_density(alpha=.3) + theme_classic() + ggtitle("nucleotide diversity northern") + theme(plot.title = element_text(hjust = 0.5)) + labs(x = expression(pi))
pi_south=ggplot(south_pi_comb, aes(x=PI, fill=region)) + geom_density(alpha=.3) + theme_classic() + ggtitle("nucleotide diversity southern") + theme(plot.title = element_text(hjust = 0.5)) + labs(x = expression(pi))

h_south=ggplot(south_popgenome_comb, aes(x=H, fill=region)) + geom_density(alpha=.3) + theme_classic() + xlim(-3,2) + ggtitle("Fay and Wu's H southern") + theme(plot.title = element_text(hjust = 0.5))
h_north=ggplot(north_popgenome_comb, aes(x=H, fill=region)) + geom_density(alpha=.3) + theme_classic() + xlim(-3,2) + ggtitle("Fay and Wu's H northern") + theme(plot.title = element_text(hjust = 0.5))

dprime_south_subset=ggplot(south_ld_subset_comb, aes(x=Dprime, fill=region)) + geom_density(alpha=.3) + theme_classic() + ggtitle("D' southern") + theme(plot.title = element_text(hjust = 0.5)) + labs(x = "D'")
dprime_north=ggplot(north_ld_comb, aes(x=Dprime, fill=region)) + geom_density(alpha=.3) + theme_classic() + ggtitle("D' northern") + theme(plot.title = element_text(hjust = 0.5)) + labs(x = "D'")

#Split into two subplots
pdf("gen_summary_stats_density.1.pdf",height=10)
print(grid.arrange(pi_south, pi_north, tajd_south,tajd_north,nrow = 4))
dev.off()

pdf("gen_summary_stats_density.2.pdf",height=10)
print(grid.arrange(h_south,h_north,dprime_south_subset,dprime_north,nrow = 4))
dev.off()
```


# Supplementary Figure 8

R code to generate Supplementary Figure 8:

```
library(ggplot2)

chrom3_start=56E6
chrom5_start=3.27E5
chrom5_end=4.36E6

input_northern=read.delim("freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.pol.counts.out",header=F)
input_southern=read.delim("freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.7_samples.pol.counts.out",header=F)

#Remove sites with zero counts
input_northern=input_northern[which(input_northern[,3]>0),]
input_southern=input_southern[which(input_southern[,3]>0),]


#Chr1
input_northern_chr1=input_northern[which(input_northern[,1]=="Scaffold19"),]
input_southern_chr1=input_southern[which(input_southern[,1]=="Scaffold19"),]

input_northern_rest=input_northern[-which(input_northern[,1]=="Scaffold19"),]
input_southern_rest=input_southern[-which(input_southern[,1]=="Scaffold19"),]

#Chr3
input_northern_chr3=input_northern_rest[which(input_northern_rest[,1]=="Scaffold61" & input_northern_rest[,2]>chrom3_start),]
input_southern_chr3=input_southern_rest[which(input_southern_rest[,1]=="Scaffold61" & input_southern_rest[,2]>chrom3_start),]

input_northern_rest=input_northern_rest[-which(input_northern_rest[,1]=="Scaffold61" & input_northern_rest[,2]>chrom3_start),]
input_southern_rest=input_southern_rest[-which(input_southern_rest[,1]=="Scaffold61" & input_southern_rest[,2]>chrom3_start),]


#Chr5
input_northern_chr5=input_northern_rest[which(input_northern_rest[,1]=="Scaffold0" & input_northern_rest[,2]<chrom5_end & input_northern_rest[,2]>chrom5_start),]
input_southern_chr5=input_southern_rest[which(input_southern_rest[,1]=="Scaffold0" & input_southern_rest[,2]<chrom5_end & input_southern_rest[,2]>chrom5_start),]

input_northern_rest=input_northern_rest[-which(input_northern_rest[,1]=="Scaffold0" & input_northern_rest[,2]<chrom5_end & input_northern_rest[,2]>chrom5_start),]
input_southern_rest=input_southern_rest[-which(input_southern_rest[,1]=="Scaffold0" & input_southern_rest[,2]<chrom5_end & input_southern_rest[,2]>chrom5_start),]


#Plot chrom 1

input_northern_chr1=cbind(rep("north",nrow(input_northern_chr1)),input_northern_chr1)
names(input_northern_chr1)[1]="pop"
input_southern_chr1=cbind(rep("south",nrow(input_southern_chr1)),input_southern_chr1)
names(input_southern_chr1)[1]="pop"
input_comb_chr1=rbind(input_northern_chr1,input_southern_chr1)
names(input_comb_chr1)[2:4]=c("scaffold","pos","count")

p1=ggplot(input_comb_chr1, aes(x=count, fill=pop)) + geom_histogram(aes(y=..density..),binwidth=0.5, position="dodge") + theme_classic() + ggtitle("Chromosome 1") + theme(plot.title = element_text(hjust = 0.5)) 


#Plot chrom 3

input_northern_chr3=cbind(rep("north",nrow(input_northern_chr3)),input_northern_chr3)
names(input_northern_chr3)[1]="pop"
input_southern_chr3=cbind(rep("south",nrow(input_southern_chr3)),input_southern_chr3)
names(input_southern_chr3)[1]="pop"
input_comb_chr3=rbind(input_northern_chr3,input_southern_chr3)
names(input_comb_chr3)[2:4]=c("scaffold","pos","count")

p2=ggplot(input_comb_chr3, aes(x=count, fill=pop)) + geom_histogram(aes(y=..density..),binwidth=0.5, position="dodge") + theme_classic() + ggtitle("Chromosome 3") + theme(plot.title = element_text(hjust = 0.5)) 


#Plot chrom 5

input_northern_chr5=cbind(rep("north",nrow(input_northern_chr5)),input_northern_chr5)
names(input_northern_chr5)[1]="pop"
input_southern_chr5=cbind(rep("south",nrow(input_southern_chr5)),input_southern_chr5)
names(input_southern_chr5)[1]="pop"
input_comb_chr5=rbind(input_northern_chr5,input_southern_chr5)
names(input_comb_chr5)[2:4]=c("scaffold","pos","count")


p3=ggplot(input_comb_chr5, aes(x=count, fill=pop)) + geom_histogram(aes(y=..density..),binwidth=0.5, position="dodge") + theme_classic() + ggtitle("Chromosome 5") + theme(plot.title = element_text(hjust = 0.5)) 


#Plot rest of genome


input_northern_rest=cbind(rep("north",nrow(input_northern_rest)),input_northern_rest)
names(input_northern_rest)[1]="pop"
input_southern_rest=cbind(rep("south",nrow(input_southern_rest)),input_southern_rest)
names(input_southern_rest)[1]="pop"
input_comb_rest=rbind(input_northern_rest,input_southern_rest)
names(input_comb_rest)[2:4]=c("scaffold","pos","count")


p4=ggplot(input_comb_rest, aes(x=count, fill=pop)) + geom_histogram(aes(y=..density..),binwidth=0.5, position="dodge") + theme_classic() + ggtitle("Rest of genome") + theme(plot.title = element_text(hjust = 0.5)) 



library(gridExtra)

pdf("allele_frequency_spectrum_ww.pdf",height=10)

print(grid.arrange(p1,p2,p3,p4,nrow = 4))

dev.off()
```

# Supplementary Figure 9

R code to plot Supplementary Figure 9:

```

ld_south=read.delim("ww_southern_hifi_bionano.filt.southern.7_samples.ld.out",header=F)
ld_north=read.delim("ww_southern_hifi_bionano.filt.northern.ld.out",header=F)

popgenome_south=read.delim("southern_with_anc_popgenome.2.out",header=F)
popgenome_north=read.delim("northern_with_anc_popgenome.2.out",header=F)

south_pi=read.delim("freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.pi_10kb.out",header=T)
north_pi=read.delim("freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.pi_10kb.out",header=T)

xpnsl=read.delim("xpnsl.10kb.out",header=T)
xpnsl=xpnsl[xpnsl$frac_scores_gt_threshold>0,]


### Scaffold0: Chromosome 5

#Sweepfinder
sweepfinder_south_scaff=read.delim("Scaffold0.southern.sweepfinder.out")
sweepfinder_north_scaff=read.delim("Scaffold0.northern.sweepfinder.out")


#XPNSL
xpnsl_scaffold=xpnsl[which(xpnsl$scaffold=="Scaffold0"),]
xpnsl_scaffold=xpnsl_scaffold[which(xpnsl_scaffold$end<5E6),]

#Fay AND Wu's H
popgenome_south_scaff=popgenome_south[which(popgenome_south[,1]=="Scaffold0"),]
popgenome_north_scaff=popgenome_north[which(popgenome_north[,1]=="Scaffold0"),]

#Pi
south_pi_scaff=south_pi[which(south_pi[,1]=="Scaffold0"),]
north_pi_scaff=north_pi[which(north_pi[,1]=="Scaffold0"),]

#LD
ld_south_scaff=ld_south[which(ld_south[,1]=="Scaffold0"),]
ld_north_scaff=ld_north[which(ld_north[,1]=="Scaffold0"),]

xstart=3E5
xend=4.4E6

pdf("chrom5.gen_summary_stats.10kb.pdf")

layout(matrix(1:5, ncol = 1), widths = 1, heights = c(0.5,0.5,0.5,0.5,0.6), respect = FALSE)
par(mar = c(0.3, 4.1, 2, 2.1))

plot(rowMeans(xpnsl_scaffold[,2:3]),xpnsl_scaffold$frac_scores_gt_threshold,ylim=c(0,1),ylab="fraction of outlier SNPs",xlim=c(xstart,xend),type="l",col="green",lwd=2,main="Chrom 5 (Scaffold0)",xaxt="n",xlab="")
points(rowMeans(xpnsl_scaffold[,2:3]),xpnsl_scaffold$frac_scores_lt_threshold,type="l",col="blue",lwd=2)
par(mar = c(0.3, 4.1, 0.5, 2.1))

plot(sweepfinder_south_scaff$location,sweepfinder_south_scaff$LR,ylim=c(0,500),xlim=c(xstart,xend),type="l",col="green",lwd=2,xaxt="n",ylab="CLR",main="")
points(sweepfinder_north_scaff$location,sweepfinder_north_scaff$LR,col="blue",type="l",lwd=2)

plot(rowMeans(popgenome_south_scaff[,2:3]),popgenome_south_scaff[,5],ylab="H",main="",xlim=c(xstart,xend),ylim=c(-6,1),type="l",lwd=2,col="green",xaxt="n",xlab="") 
points(rowMeans(popgenome_north_scaff[,2:3]),popgenome_north_scaff[,5],type="l",lwd=2,col="blue") 

plot(rowMeans(south_pi_scaff[,2:3]),south_pi_scaff[,5],xlab="position",ylab=expression(pi),main="",xlim=c(xstart,xend),ylim=c(0,0.017),type="l",lwd=2,col="green",xaxt="n") 
points(rowMeans(north_pi_scaff[,2:3]),north_pi_scaff[,5],type="l",lwd=2,col="blue")

par(mar = c(4.1, 4.1, 0.5, 2.1))
plot(rowMeans(ld_south_scaff[,2:3]),ld_south_scaff[,5],xlab="position",ylab="D'",xlim=c(xstart,xend),ylim=c(0,1),type="l",lwd=2,main="",col="green") 
points(rowMeans(ld_north_scaff[,2:3]),ld_north_scaff[,5],xlab="position",type="l",lwd=2,col="blue")
dev.off()


### Scaffold19 - chrom 1

#Sweepfinder data
sweepfinder_south_scaff=read.delim("Scaffold19.southern.sweepfinder.out")
sweepfinder_north_scaff=read.delim("Scaffold19.northern.sweepfinder.out")

xpnsl_scaffold=xpnsl[which(xpnsl$scaffold=="Scaffold19"),]

#Fay and Wu's H
popgenome_south_scaff=popgenome_south[which(popgenome_south[,1]=="Scaffold19"),]
popgenome_north_scaff=popgenome_north[which(popgenome_north[,1]=="Scaffold19"),]

#Pi
south_pi_scaff=south_pi[which(south_pi[,1]=="Scaffold19"),]
north_pi_scaff=north_pi[which(north_pi[,1]=="Scaffold19"),]

#LD
ld_south_scaff=ld_south[which(ld_south[,1]=="Scaffold19"),]
ld_north_scaff=ld_north[which(ld_north[,1]=="Scaffold19"),]

pdf("chrom1.gen_summary_stats.10kb.pdf")

layout(matrix(1:5, ncol = 1), widths = 1, heights = c(0.5,0.5,0.5,0.5,0.6), respect = FALSE)

par(mar = c(0.3, 4.1, 2, 2.1))
plot(rowMeans(xpnsl_scaffold[,2:3]),xpnsl_scaffold$frac_scores_gt_threshold,ylim=c(0,1),xlab="",ylab="fraction outlier SNPs",type="l",col="green",lwd=2,main="Chrom 1 (Scaffold19)",xaxt="n")
points(rowMeans(xpnsl_scaffold[,2:3]),xpnsl_scaffold$frac_scores_lt_threshold,type="l",col="blue",lwd=2)

par(mar = c(0.3, 4.1, 0.5, 2.1))
plot(sweepfinder_south_scaff$location,sweepfinder_south_scaff$LR,ylim=c(0,500),type="l",col="green",lwd=2,xaxt="n",ylab="CLR",main="")
points(sweepfinder_north_scaff$location,sweepfinder_north_scaff$LR,col="blue",type="l",lwd=2)

plot(rowMeans(popgenome_south_scaff[,2:3]),popgenome_south_scaff[,5],xlab="",ylab="H",main="",ylim=c(-6,1),type="l",lwd=2,col="green",xaxt="n") 
points(rowMeans(popgenome_north_scaff[,2:3]),popgenome_north_scaff[,5],type="l",lwd=2,col="blue") 

plot(rowMeans(south_pi_scaff[,2:3]),south_pi_scaff[,5],xlab="position",ylab=expression(pi),main="",type="l",lwd=2,col="green",xaxt="n",ylim=c(0,0.017)) 
points(rowMeans(north_pi_scaff[,2:3]),north_pi_scaff[,5],type="l",lwd=2,col="blue") 

par(mar = c(4.1, 4.1, 0.5, 2.1))
plot(rowMeans(ld_south_scaff[,2:3]),ld_south_scaff[,5],xlab="position",ylab="D'",ylim=c(0,1),type="l",lwd=2,main="",col="green") 
points(rowMeans(ld_north_scaff[,2:3]),ld_north_scaff[,5],xlab="position",type="l",lwd=2,col="blue")
dev.off()



### Scaffold61 - chrom 3 

#Sweepfinder
sweepfinder_south_scaff=read.delim("Scaffold61.southern.sweepfinder.out")
sweepfinder_north_scaff=read.delim("Scaffold61.northern.sweepfinder.out")


#Xpnsl
xpnsl_scaffold=xpnsl[which(xpnsl$scaffold=="Scaffold61"),]


#Fay and Wu's H
popgenome_south_scaff=popgenome_south[which(popgenome_south[,1]=="Scaffold61"),]
popgenome_north_scaff=popgenome_north[which(popgenome_north[,1]=="Scaffold61"),]


#Pi
south_pi_scaff=south_pi[which(south_pi[,1]=="Scaffold61"),]
north_pi_scaff=north_pi[which(north_pi[,1]=="Scaffold61"),]

#LD
ld_south_scaff=ld_south[which(ld_south[,1]=="Scaffold61"),]
ld_north_scaff=ld_north[which(ld_north[,1]=="Scaffold61"),]


xstart=56E6
xend=69.5E6


pdf("chrom3.gen_summary_stats.10kb.pdf")

layout(matrix(1:5, ncol = 1), widths = 1, heights = c(0.5,0.5,0.5,0.5,0.6), respect = FALSE)

par(mar = c(0.3, 4.1, 2, 2.1))
plot(rowMeans(xpnsl_scaffold[,2:3]),xpnsl_scaffold$frac_scores_gt_threshold,ylim=c(0,1),xlab="",ylab="fraction outlier SNPs",type="l",col="green",lwd=2,xlim=c(56E6,70E6),main="Chrom 3 (Scaffold61)",xaxt="n")
points(rowMeans(xpnsl_scaffold[,2:3]),xpnsl_scaffold$frac_scores_lt_threshold,ylim=c(0,0.5),type="l",col="blue",lwd=2,xlim=c(xstart,xend))

par(mar = c(0.3, 4.1, 0.5, 2.1))
plot(sweepfinder_south_scaff$location,sweepfinder_south_scaff$LR,ylim=c(0,500),xlim=c(xstart,xend),type="l",col="green",lwd=2,xaxt="n",ylab="CLR",main="")
points(sweepfinder_north_scaff$location,sweepfinder_north_scaff$LR,col="blue",type="l",lwd=2)

plot(rowMeans(popgenome_south_scaff[,2:3]),popgenome_south_scaff[,5],xlab="",ylab="H",main="",ylim=c(-6,1),xlim=c(xstart,xend),type="l",lwd=2,col="green",xaxt="n") 
points(rowMeans(popgenome_north_scaff[,2:3]),popgenome_north_scaff[,5],type="l",lwd=2,col="blue") 

plot(rowMeans(south_pi_scaff[,2:3]),south_pi_scaff[,5],xlab="",ylab=expression(pi),main="",type="l",lwd=2,col="green",xlim=c(xstart,xend),ylim=c(0,0.017),xaxt="n") 
points(rowMeans(north_pi_scaff[,2:3]),north_pi_scaff[,5],type="l",lwd=2,col="blue") 

par(mar = c(4.1, 4.1, 0.5, 2.1))
plot(rowMeans(ld_south_scaff[,2:3]),ld_south_scaff[,5],xlab="position",ylab="D'",ylim=c(0,1),type="l",lwd=2,main="",col="green",xlim=c(xstart,xend)) 
points(rowMeans(ld_north_scaff[,2:3]),ld_north_scaff[,5],xlab="position",type="l",lwd=2,col="blue")

dev.off()
```

# Supplementary Figure 10

R code to generate Supplementary Figure 10 - selection signatures in the SPON1 gene


```
#Read gene data, which are simply SPON1 exon and CDS intervals taken from the gene annotation file
gene_data=read.delim("/proj/snic2020-16-64/max/hifi_data/assemblies/annotation_data/SPON1.intervals.txt",header=F)


#Read diversity data
south_pi=read.delim("freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.southern.pi_10kb.out",header=T)
north_pi=read.delim("freebayes_reseq_southern_hifi_bionano.filt3.decomposed.no_repeat.new.phased.bi_allelic_snps.northern.pi_10kb.out",header=T)


#Read XP-nsl data
xpnsl=read.delim("xpnsl.10kb.out",header=T)
xpnsl=xpnsl[xpnsl$frac_scores_gt_threshold>=0,]

xpnsl_scaff=xpnsl[which(xpnsl[,1]=="Scaffold0"),]


#Read Sweepfinder2 data
sweepfinder_south=read.delim("Scaffold0.southern.sweepfinder.out",header=T)
sweepfinder_north=read.delim("Scaffold0.northern.sweepfinder.out",header=T)


pdf("Figure3.pdf")

par(mar = c(0.3, 4.1, 0.5, 2.1))

layout(matrix(1:4, ncol = 1), widths = 1, heights = c(0.5,0.5,0.5,0.5), respect = FALSE)

plot(rowMeans(xpnsl_scaff[,2:3]),xpnsl_scaff$frac_scores_gt_threshold,ylim=c(0,1),ylab="fraction of outlier SNPs",type="l",col="green",lwd=2,xlim=c(5.25e5,7.5E5),xaxt="n",xlab="")
points(rowMeans(xpnsl_scaff[,2:3]),xpnsl_scaff$frac_scores_lt_threshold,type="l",col="blue",lwd=2)

plot(sweepfinder_south$location,sweepfinder_south$LR,col="green",ylab="CLR",xlim=c(5.25e5,7.5E5),type="l",main="",xaxt="n",ylim=c(0,500),lwd=2)
points(sweepfinder_north$location,sweepfinder_north$LR,col="blue",type="l",lwd=2)

plot(rowMeans(south_pi_scaff[,2:3]),south_pi_scaff$PI,ylim=c(0,0.01),ylab=expression(pi),type="l",col="green",lwd=2,xlim=c(5.25e5,7.5E5),xaxt="n",xlab="")
points(rowMeans(north_pi_scaff[,2:3]),north_pi_scaff$PI,xlab="position",type="l",col="blue",lwd=2)

par(mar = c(4.1, 4.1, 0.3, 2.1))

plot(NULL,yaxt="n", xlab="position",xlim=c(5.25e5,7.5E5),ylim=c(-0.6,-0.2))

exondata=gene_data[which(gene_data[,2]=="exon"),]
cdsdata=gene_data[which(gene_data[,2]=="CDS"),]
gene_color="red"

segments(x0=min(exondata[,3]),x1=max(exondata[,4]),y0=-0.325,y1=-0.325,lwd=2,col=gene_color,border=gene_color)

upper_exon=-0.30
lower_exon=-0.35

for(i in c(1:nrow(exondata))){

start=exondata[i,3]
end=exondata[i,4]

polygon(x=c(start,end,end,start),y=c(lower_exon,lower_exon,upper_exon,upper_exon),col=gene_color,border=gene_color)

}


upper_cds=-0.275
lower_cds=-0.375

for(i in c(1:nrow(cdsdata))){

start=cdsdata[i,3]
end=cdsdata[i,4]

polygon(x=c(start,end,end,start),y=c(lower_cds,lower_cds,upper_cds,upper_cds),col=gene_color,border=gene_color)

}

text(y=-0.425,x=mean(c(min(exondata[,3]),max(exondata[,4])))+500,"SPON1")


dev.off()
```

