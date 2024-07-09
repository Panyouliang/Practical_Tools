#https://www.genome.jp/kegg-bin/get_htext?ko00001

if(!file.exists('kegg_info.RData')){
  
  library(jsonlite)
  library(purrr)
  library(RCurl)
  library(stringr)
  library(tibble)
  
  update_kegg <- function(json = "ko00001.json",file=NULL) {
    #pathway2name <- tibble(Pathway = character(), Name = character())
    pathway2name <- tibble(Lev1= character(),Lev2=character(),Pathway = character(), Name = character())
    ko2pathway <- tibble(Ko = character(), Pathway = character())
    
    kegg <- fromJSON(json)
    
    for (a in seq_along(kegg[["children"]][["children"]])) {

      A <- kegg[["children"]][["name"]][[a]]
      for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
        B <- kegg[["children"]][["children"]][[a]][["name"]][[b]] 
        
        for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
          pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
          
          pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
          pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5} ", "")
          #pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Name = pathway_name))
          pathway2name <- rbind(pathway2name, tibble(Lev1 = A, Lev2=B ,Pathway = pathway_id, Name = pathway_name))
          
          kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
          
          kos <- str_match(kos_info, "K[0-9]*")[,1]
          
          ko2pathway <- rbind(ko2pathway, tibble(Ko = kos, Pathway = rep(pathway_id, length(kos))))
        }
      }
    }
    
    save(pathway2name, ko2pathway, file = file)
  }
  
  update_kegg(json = "ko00001.json",file="kegg_info.RData")
  
}

#load("kegg_info.RData")
library(rjson)
my.JSON <- fromJSON(file="ko00002.json")
x=my.JSON[2]
df <- as.data.frame(do.call(rbind, x))
n=data.frame(matrix(ncol = 4))
colnames(n)=c("Type","Lev1","Lev2","Lev3")
for (i in seq(1:2)) {
  #i=1
  a=df[1,i]
  a.df= as.data.frame(do.call(rbind, a))
  type=unlist(a.df[1,1])
  b=a.df[1,2]
  b.df=as.data.frame(do.call(rbind, b))
  for (j in seq(1:ncol(b.df))) {
    #j=1
    c=b.df[1,j]
    c.df=as.data.frame(do.call(rbind, c))
    Lev1=unlist(c.df[1,1])
    d=as.data.frame(do.call(rbind, c.df[1,2]))
    for (x in seq(1:ncol(d))) {
      #x=1
      e=d[1,x]
      e.df=as.data.frame(do.call(rbind, e))
      Lev2=unlist(e.df[1,1])
      f=as.data.frame(do.call(rbind, e.df[1,2]))
      for (y in seq(1:ncol(f))) {
        #y=1
        g=f[1,y]
        g.df= as.data.frame(do.call(rbind, g))
        Lev3=unlist(g.df[1,1])
        t=data.frame(Type=type,Lev1=Lev1,Lev2=Lev2,Lev3=Lev3)
        n=rbind(n,t)
      }
    }
  }
}
n=n[-1,]

library(stringr)
library(tidyr)
nl=separate(n,Lev3,into = c("Module","Pathway"),sep = "\\[PATH:")
nl$Pathway=sub("\\]","",nl$Pathway)
nlw=nl %>% separate_rows(Pathway, sep = ' ', convert = F)
nlw$Pathway=sub("map","ko",nlw$Pathway)
load("kegg_info.RData")
keg_lvs=merge(pathway2name,nlw[,c(4,5)],by="Pathway",all = T)
a=keg_lvs[which(is.na(keg_lvs$Lev1)),]

save(keg_lvs,ko2pathway,file="KEGG_levs.pathway.module.Rdata")
