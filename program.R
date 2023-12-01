library(maftools)
library(stringr)
sarc <- read.maf(maf = "TCGA.SARC.mutect.cc207fe8-ee0a-4b65-82cb-c8197d264126.DR-10.0.somatic.maf")
load("prot_db.RData")

Gene_Name="TP53"
Ref_seq="NM_001126112"


#Draw Entier Gene. it shows in a gray rectangle

index=which(prot_db[,2]==Ref_seq)         #Find The Ref_seq In Input Data
Gene_Index=which(prot_db[,1]==Gene_Name)  #Find The Max Gene Lenght 
Gene_lenght=max(prot_db[Gene_Index,4])

HeightOfRect=1                            #Customize the Plot Coordinate
min_ylim=-3
max_ylim=15

par(mar=c(7,5,3,1),xpd=TRUE)
plot(0:Gene_lenght,0:Gene_lenght,type="n",axes=T,ann=T,ylim = c(min_ylim,max_ylim),xlab = "",ylab = "#Mutations")
rect(0,(min_ylim+.5),Gene_lenght,-.5,col = "gray")

  
# Find The Mutations

  # Customize Plot
  All_Mutations_Types=c("Missense_Mutation","Frame_Shift_Del","Splice_Site","Nonsense_Mutation","Nonstop_Mutation",
                        "Frame_Shift_Ins","Translation_Start_Site","In_Frame_Del","In_Frame_Ins")
  Mutations_Colors=c("#2ec4b6","#43aa8b","#90be6d","#f9c74f","#f8961e","#f3722c","#f94144","#f19c79","#577590")
  j=1
  legend_label=0
  legend_Color=0
  
  # Searching For Mutations
  index=which(sarc@data[["Hugo_Symbol"]]==Gene_Name)                            #1: Finding Mutation index in the input gene
  Mutations=sarc@data[["Variant_Classification"]][index]                        #2: Finding Type Of The Mutations
  Mutations_Type = str_split(Mutations,"\n")
  Somatic_Mutation_Rate=ceiling((length(Mutations_Type)/Gene_lenght)*10000)     #3: Calculating Somatic Mutation Rate
  Somatic_Mutation_Rate=Somatic_Mutation_Rate/100
  for(i in 1:9)                                                                 #4: Finding Mutations Locations And Add Them In To The Plot
  {
      Grouped=which(Mutations_Type==All_Mutations_Types[i])
    if(length(Grouped)>0)
    {
      legend_label[j]=All_Mutations_Types[i]
      
      Mutations_Locations = sarc@data[["HGVSp_Short"]][index[Grouped]]
      Mutations_Locations_X=as.integer(str_extract(Mutations_Locations, "\\d+")) 
      Mutations_Locations_X=sort(Mutations_Locations_X)
      
      No_Of_Mutation_In_Gene=length(Mutations_Locations_X)
      
      Unique=unique(Mutations_Locations_X)
      temp=length(unique(Mutations_Locations_X))
      
      while(temp>0)
      {
        No_Of_Mutation_In_ALocation=length(which(Mutations_Locations_X==Unique[temp]))
        points(Unique[temp],No_Of_Mutation_In_ALocation,pch=19 , col=Mutations_Colors[i],cex=1.5)
        segments(Unique[temp],-.5,Unique[temp],No_Of_Mutation_In_ALocation,col = Mutations_Colors[i])
        temp=temp-1
        
      }
      
      legend_Color[j]=Mutations_Colors[i]
      j=j+1
    }
  }


  
  # Finding Domains
  
  index=which(prot_db[,2]==Ref_seq)       #Search Domains In The prot_db
  
  Domains_Start=prot_db[index,5]          #Finding Domains Start Locations
  Domains_Start=Domains_Start[[1]]
  
  Domains_End=prot_db[index,6]            #Finding Domains End Locations
  Domains_End=Domains_End[[1]]
  
  Number_Of_Domains=length(Domains_Start) #Calculating Number Of Domains In A Gene
  
  
  # Plot Domains And Add Them To Legend

  PlotColors=c("#577590","#43aa8b","#90be6d","#f9c74f","#f8961e","#f3722c","#f94144");
  
  Domain_Label=prot_db[index,8]
  Domain_Label=Domain_Label[[1]]
  Domain_Label_Location=floor((Domains_Start+Domains_End)/2)-5
  
  Domains_Colors=PlotColors[1:(Number_Of_Domains+1)]
  rect(Domains_Start,min_ylim,Domains_End,0,col=Domains_Colors)
  
  text(Domain_Label_Location,y=min_ylim/2,Domain_Label,adj = 0,cex=.7)
  legend("bottomleft",pch=20,inset=c(0,-.5),legend=legend_label,fill=legend_Color,cex=0.55,box.lty=0,y.intersp=0.5,x.intersp =.5,ncol=5, nrow(1) )
  legend("topleft",inset = c(-0.04,-.2),legend=(paste(Gene_Name," :[Somatic_Mutation_Rate: ",Somatic_Mutation_Rate,"%]","\n",Ref_seq, collapse=""))
         ,cex=1,box.lty = 0,y.intersp=.1,x.intersp =0)
  

  