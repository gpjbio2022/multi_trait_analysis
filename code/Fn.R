##FN00

library(tidyverse)
library(skimr)

# read in data
Jf<-read_csv("~/Desktop/精分.csv")

Zb<-read_csv("~/Desktop/自闭.csv")

Zy<-read_csv("~/Desktop/躁郁.csv")

# join three tables
J3Dt<-full_join(x=Jf,y=Zy,by="Gene")%>%
  full_join(y=Zb,by="Gene")

# EDA Jf
Jf

skim(Jf)

##rename

Jfr <- Jf %>%
  rename(CasePTV=`Case PTV`,
         ControlPTV=`Control PTV`,
         CaseMissense_B3=`Case Missense (MPC ≥ 3)`,
         ControlMissense_B3=`Control Missense (MPC ≥ 3)`,
         CaseMissense_B23=`Case Missense (2 ≤ MPC < 3)`,
         ControlMissense_B23=`Control Missense (2 ≤ MPC < 3)`,
         DenovoPTV=`De novo PTV`,
         DenovoMissense_B3=`De Novo Missense (MPC ≥ 3)`,
         DenovoMissense_B23=`De Novo Missense (2 ≤ MPC < 3)`,
         Pm=`P meta`,
         Qm=`Q meta`,
         OR_I=`OR (Class I)`,
         OR_II=`OR (Class II)`)

Zbr <- Zb %>%
  rename(CasePTV=`Case PTV`,
         ControlPTV=`Control PTV`,
         CaseMissense_B3=`Case Missense (MPC ≥ 3)`,
         ControlMissense_B3=`Control Missense (MPC ≥ 3)`,
         CaseMissense_B23=`Case Missense (2 ≤ MPC < 3)`,
         ControlMissense_B23=`Control Missense (2 ≤ MPC < 3)`,
         DenovoPTV=`De novo PTV`,
         DenovoMissense_B3=`De Novo Missense (MPC ≥ 3)`,
         DenovoMissense_B23=`De Novo Missense (2 ≤ MPC < 3)`,
         Pm=`P meta`,
         Qm=`Q meta`,
         OR_I=`OR (Class I)`,
         OR_II=`OR (Class II)`)


Zyr <- Zy %>%
  rename(CasePTV=`Case PTV`,
         ControlPTV=`Control PTV`,
         CaseMissense_B3=`Case Missense (MPC ≥ 3)`,
         ControlMissense_B3=`Control Missense (MPC ≥ 3)`,
         CaseMissense_B23=`Case Missense (2 ≤ MPC < 3)`,
         ControlMissense_B23=`Control Missense (2 ≤ MPC < 3)`,
         DenovoPTV=`De novo PTV`,
         DenovoMissense_B3=`De Novo Missense (MPC ≥ 3)`,
         DenovoMissense_B23=`De Novo Missense (2 ≤ MPC < 3)`,
         Pm=`P meta`,
         Qm=`Q meta`,
         OR_I=`OR (Class I)`,
         OR_II=`OR (Class II)`)

ggplot(data=Jf_rename)+geom_bar(mapping=aes(x=Case_PTV))

Jf %>% 
  sample_frac(0.2) %>%
  ggplot()+
  geom_point(mapping=aes(x=ControlMissense_B3,y=CaseMissense_B3),position="jitter",alpha=0.8,color="blue",show.legend=TRUE)+
  geom_point(mapping=aes(x=ControlMissense_B23,y=CaseMissense_B23),position="jitter",alpha=0.2,color="red",show.legend=TRUE)+
  geom_abline(color="green")

Jf_long <- Jf %>%qw
pivot_longer(3:8, names_to = "type", values_to = "count")

# Jf_long %>%
#   filter(type == "CasePTV" | type == "ControlPTV") %>%
#   ggplot() +
#    geom_histogram(aes(type, group="type"))
#
# geom_point(
#   if(type=="CasePTV"){
#     color="red"
#   }
#   if(type=="Case"){...}
# )

?case_when
Jf_long %>%
  mutate(color = case_when(
    type == "CasePTV" ~ "red",
    type == "ControlPTV" ~ "blue",
    TRUE ~ "white"
  )) %>% view()



#ALL JOIN +后缀

all_join()

?full_join

# fisher exact test
?fisher.test

TeaTasting <-
  matrix(c(3, 1, 1, 3),
         nrow = 2,
         dimnames = list(Guess = c("Milk", "Tea"),
                         Truth = c("Milk", "Tea")))

TeaTasting

fisher.test(TeaTasting, alternative = "two.sided")

EX1<-matrix(c(16,0,13917,14422),
            nrow=2,
            dimnames=list(Grouptype=c("Case", "Control"),
                          Datatype=c("PTV", "NonePTV")))
EX1.fisher <- fisher.test(EX1, alternative = "two.sided")

str(EX1.fisher)

EX1.fisher$p.value
EX1.fisher$estimate

JfTemp <- Jf %>%
  replace_na(list(CasePTV = 0, ControlPTV = 0))
skim(JfTemp)

Jf.test <- Jf %>%
  mutate(OPTV = fisher.test(matrix(c(CasePTV, ControlPTV,(13933-CasePTV),(14422-ControlPTV)), nrow = 2, byrow = T))$estimate,
         PPTV = fisher.test(matrix(c(CasePTV, ControlPTV,(13933-CasePTV),(14422-ControlPTV)), nrow = 2, byrow = T))$p.value)

INFrplc<-function(x) ifelse(is.infinite(x),114514,x)

lJf<-filter(Jf,!is.na(CasePTV)&&!is.na(ControlPTV))

Jf.test <- Jf %>%
  mutate(OPTV=(unname(INFrplc(fisher.test(matrix(c(CasePTV, ControlPTV,(13933-CasePTV),(14422-ControlPTV)), nrow = 2, byrow = T))$estimate))),
         PPTV=(unname(INFrplc(fisher.test(matrix(c(CasePTV, ControlPTV,(13933-CasePTV),(14422-ControlPTV)), nrow = 2, byrow = T))$p.value))))

Jf.test <- Jf

mutate(Jf.test,OPTV = as.numeric(INFrplc(fisher.test(matrix(c(CasePTV, ControlPTV,(13933-CasePTV),(14422-ControlPTV)), nrow = 2, byrow = T))$estimate)))

asd<-select(Jf,CasePTV)

mutate(asd,(unname(INFrplc(fisher.test(matrix(c(CasePTV, CasePTV,(13933),(14422)), nrow = 2, byrow = T))$p.value))))
mutate(asd,rm=t.test(c(CasePTV,CasePTV,13917-CasePTV,14422-CasePTV))$p.value)


Jf.test <- Jf %>%
  mutate(OPTV = t.test(c(CasePTV, ControlPTV,(13933-CasePTV),(14422-ControlPTV)))$estimate,
         PPTV = t.test(c(CasePTV, ControlPTV,(13933-CasePTV),(14422-ControlPTV)))$p.value)

Jf.test <- Jf %>%
  filter(!is.na(CasePTV))%>%
  filter(!is.na(ControlPTV))%>%
  rowwise()%>%
  mutate(OPTV=(CasePTV/(13833-CasePTV))/(ControlPTV/(13833-ControlPTV)),
         PPTV=t.test(c(CasePTV, ControlPTV,(13933-CasePTV),(14422-ControlPTV)))$p.value)

###functionlize

Jf.test <- Jf %>%
  filter(!is.na(CasePTV))%>%
  filter(!is.na(ControlPTV))%>%
  rowwise()%>%
  mutate(OPTV=(CasePTV/(13833-CasePTV))/(ControlPTV/(13833-ControlPTV)),
         PPTV=(unname(INFrplc(fisher.test(matrix(c(CasePTV, ControlPTV,(13933-CasePTV),(14422-ControlPTV)), nrow = 2, byrow = T))$p.value))))





fisher.or <- function(a,b,m,n){
  data <- matrix(c(a,b,m-a,n-b), byrow = T, ncol = 2)
  fisher.test(data)$estimate[[1]]
}

# a function to extract Pvalue from fisher.test
fisher.pval <- function(a,b,m,n){
  data <- matrix(c(a,b,m-a,n-b), byrow = T, ncol = 2)
  fisher.test(data)$p.value
}

# first replace na for missing count number, then run fisher.test
bipolar_new <- bipolar %>%
  replace_na(list(`PTV Case Count` = 0,
                  `PTV Control Count` = 0)) %>%
  mutate(OR_PTV = pmap(list(`PTV Case Count`, `PTV Control Count`, Cases, Controls), fisher.or),
         Pval_PTV = pmap(list(`PTV Case Count`, `PTV Control Count`, Cases, Controls), fisher.pval))
##join


J3Dt<-full_join(x=Jf,y=Zy,by="Gene")%>%
  full_join(y=Zb,by="Gene")










#TYPE



fisher.or <- function(a,b,m,n){
  data <- matrix(c(a,b,m-a,n-b), byrow = T, ncol = 2)
  fisher.test(data)$estimate[[1]]
}

fisher.pval <- function(a,b,m,n){
  data <- matrix(c(a,b,m-a,n-b), byrow = T, ncol = 2)
  fisher.test(data)$p.value[[1]]
}

# first replace na for missing count number, then run fisher.test

nJF <- Jf %>%
  replace_na(list(CasePTV_Jf = 0,
                  ControlPTV_Jf = 0)) %>%
  mutate(OR_PTV = pmap(list(CasePTV_Jf, ControlPTV_Jf, 24248, 97322), fisher.or),
         Pval_PTV = pmap(list(CasePTV_Jf, ControlPTV_Jf, 24248, 97322), fisher.pval))


QPcal<-function(Jf,CasePTV_Jf,ControlPTV_Jf,maxa,maxb,adda,addb){
  Jf<-Jf %>%
    replace_na(list(CasePTV_Jf = 0,
                    ControlPTV_Jf = 0)) %>%
    mutate(!!adda:=pmap(list(CasePTV_Jf, ControlPTV_Jf, maxa, maxb),fisher.or),
           !!addb:=pmap(list(CasePTV_Jf, ControlPTV_Jf, maxa, maxb),fisher.pval))
  return(Jf)
}


QPcal<-function(Target,ExpCase,ExpControl,Case,Control,orN,pvalN){
  returning<-Target %>%
    #replace_na(list((!!sym(ExpCase)) = 0,
    #                (!!sym(ExpControl)) = 0)) ##%>%
    # mutate(!!orN:=pmap(list((!!sym(ExpCase)),(!!sym(ExpControl)),Case,Control),fisher.or),
    #        !!pvalN:=pmap(list((!!sym(ExpCase)),(!!sym(ExpControl)),Case,Control),fisher.pval))
    return(returning)
}

ns<-QPcal(nJF,CasePTV_Jf,ControlPTV_Jf,24248,97322,"OR_PTVs","Pval_PTVs")

#自我赋值内容




QPcal<-function(Target,ExpCase,ExpControl,Case,Control,orN,pvalN){
  returning<-Target %>%
    replace_na(
      list(
        (!!sym(ExpCase)) := 0,
        (!!sym(ExpControl)) := 0
      )
    )%>%
    mutate(!!orN:=pmap(list((!!sym(ExpCase)),(!!sym(ExpControl)),Case,Control),fisher.or),
           !!pvalN:=pmap(list((!!sym(ExpCase)),(!!sym(ExpControl)),Case,Control),fisher.pval))
  return(returning)
}




###

# Zb<-mutate(Zb,CasePTV_Zb=CasePTVSWE_Zb+CasePTVDEN_Zb)
# 
# Zb<-mutate(Zb,ControlPTV_Zb=ControlPTVSWE_Zb+ControlPTVDEN_Zb)
# 
# Zb<-Zb %>%
#   replace_na(list(CaseDenovoPTV_Zb = 0,
#                   ControlDenovoPTV_Zb = 0)) %>%
#   mutate(ODenovoPTV_Zb=pmap(list(CaseDenovoPTV_Zb, ControlDenovoPTV_Zb, 11257, 18501),fisher.or),
#          PDenovoPTV_Zb=pmap(list(CaseDenovoPTV_Zb, ControlDenovoPTV_Zb, 11257, 18501),fisher.pval))
# 
# Zb<-Zb %>%
#   replace_na(list(CaseDenovoMisA_Zb = 0,
#                   ControlDenovoMisA_Zb = 0)) %>%
#   mutate(ODenovoMisA_Zb=pmap(list(CaseDenovoMisA_Zb, ControlDenovoMisA_Zb, 11257, 18501),fisher.or),
#          PDenovoMisA_Zb=pmap(list(CaseDenovoMisA_Zb, ControlDenovoMisA_Zb, 11257, 18501),fisher.pval))
# 
# Zb<-Zb %>%
#   replace_na(list(CaseDenovoMisB_Zb = 0,
#                   ControlDenovoMisB_Zb = 0)) %>%
#   mutate(ODenovoMisB_Zb=pmap(list(CaseDenovoMisB_Zb, ControlDenovoMisB_Zb, 11257, 18501),fisher.or),
#          PDenovoMisB_Zb=pmap(list(CaseDenovoMisB_Zb, ControlDenovoMisB_Zb, 11257, 18501),fisher.pval))
# 
# Jf<-Jf %>%
#   replace_na(list(CaseMissense_B23_Jf = 0,
#                   ControlMissense_B23_Jf = 0)) %>%
#   mutate(OMissense_B23_Jf=pmap(list(CaseMissense_B23_Jf, ControlMissense_B23_Jf, 24248, 97322),fisher.or),
#          PMissense_B23_Jf=pmap(list(CaseMissense_B23_Jf, ControlMissense_B23_Jf, 24248, 97322),fisher.pval))
# 
# Jf<-Jf %>%
#   replace_na(list(CaseMissense_B3_Jf = 0,
#                   ControlMissense_B3_Jf = 0)) %>%
#   mutate(OMissense_B3_Jf=pmap(list(CaseMissense_B3_Jf, ControlMissense_B3_Jf, 24248, 97322),fisher.or),
#          PMissense_B3_Jf=pmap(list(CaseMissense_B3_Jf, ControlMissense_B3_Jf, 24248, 97322),fisher.pval))
# 
# Zb<-mutate(Zb,CasePTV_Zb=CasePTVSWE_Zb+CasePTVDEN_Zb)
# 
# Zb<-mutate(Zb,ControlPTV_Zb=ControlPTVSWE_Zb+ControlPTVDEN_Zb)
# 
# Zb<-Zb %>%
#   replace_na(list(CaseDenovoPTV_Zb = 0,
#                   ControlDenovoPTV_Zb = 0)) %>%
#   mutate(ODenovoPTV_Zb=as.numeric(pmap(list(CaseDenovoPTV_Zb, ControlDenovoPTV_Zb, 11257, 18501),fisher.or)),
#          PDenovoPTV_Zb=as.numeric(pmap(list(CaseDenovoPTV_Zb, ControlDenovoPTV_Zb, 11257, 18501),fisher.pval)))
# 
# Zb<-Zb %>%
#   replace_na(list(CaseDenovoMisA_Zb = 0,
#                   ControlDenovoMisA_Zb = 0)) %>%
#   mutate(ODenovoMisA_Zb=as.numeric(pmap(list(CaseDenovoMisA_Zb, ControlDenovoMisA_Zb, 11257, 18501),fisher.or)),
#          PDenovoMisA_Zb=as.numeric(pmap(list(CaseDenovoMisA_Zb, ControlDenovoMisA_Zb, 11257, 18501),fisher.pval)))
# 
# Zb<-Zb %>%
#   replace_na(list(CaseDenovoMisB_Zb = 0,
#                   ControlDenovoMisB_Zb = 0)) %>%
#   mutate(ODenovoMisB_Zb=as.numeric(pmap(list(CaseDenovoMisB_Zb, ControlDenovoMisB_Zb, 11257, 18501),fisher.or)),
#          PDenovoMisB_Zb=as.numeric(pmap(list(CaseDenovoMisB_Zb, ControlDenovoMisB_Zb, 11257, 18501),fisher.pval)))
# 
# Jf<-Jf %>%
#   replace_na(list(CaseMissense_B23_Jf = 0,
#                   ControlMissense_B23_Jf = 0)) %>%
#   mutate(OMissense_B23_Jf=as.numeric(pmap(list(CaseMissense_B23_Jf, ControlMissense_B23_Jf, 24248, 97322),fisher.or)),
#          PMissense_B23_Jf=as.numeric(pmap(list(CaseMissense_B23_Jf, ControlMissense_B23_Jf, 24248, 97322),fisher.pval)))
# 
# Jf<-Jf %>%
#   replace_na(list(CaseMissense_B3_Jf = 0,
#                   ControlMissense_B3_Jf = 0)) %>%
#   mutate(OMissense_B3_Jf=as.numeric(pmap(list(CaseMissense_B3_Jf, ControlMissense_B3_Jf, 24248, 97322),fisher.or)),
#          PMissense_B3_Jf=as.numeric(pmap(list(CaseMissense_B3_Jf, ControlMissense_B3_Jf, 24248, 97322),fisher.pval)))


Zb<-mutate(Zb,CasePTV_Zb=CasePTVSWE_Zb+CasePTVDEN_Zb)

Zb<-mutate(Zb,ControlPTV_Zb=ControlPTVSWE_Zb+ControlPTVDEN_Zb)

Jf<-Jf %>%
  replace_na(list(CaseMissense_B23_Jf = 0,
                  ControlMissense_B23_Jf = 0)) %>%
  mutate(OMissense_B23_Jf=as.numeric(pmap(list(CaseMissense_B23_Jf, ControlMissense_B23_Jf, 24248, 97322),fisher.or)),
         PMissense_B23_Jf=as.numeric(pmap(list(CaseMissense_B23_Jf, ControlMissense_B23_Jf, 24248, 97322),fisher.pval)))

Jf<-Jf %>%
  replace_na(list(CaseMissense_B3_Jf = 0,
                  ControlMissense_B3_Jf = 0)) %>%
  mutate(OMissense_B3_Jf=as.numeric(pmap(list(CaseMissense_B3_Jf, ControlMissense_B3_Jf, 24248, 97322),fisher.or)),
         PMissense_B3_Jf=as.numeric(pmap(list(CaseMissense_B3_Jf, ControlMissense_B3_Jf, 24248, 97322),fisher.pval)))

Zy<-Zy %>%
  replace_na(list(CasePTV_Zy = 0,
                  ControlPTV_Zy = 0)) %>%
  mutate(tOPTV_Zy=as.numeric(pmap(list(CasePTV_Zy, ControlPTV_Zy, 13922, 14422),fisher.or)),
         tPPTV_Zy=as.numeric(pmap(list(CasePTV_Zy, ControlPTV_Zy, 13922, 14422),fisher.pval)))

Zb<-Zb %>%
  replace_na(list(CasePTV_Zb = 0,
                  ControlPTV_Zb = 0)) %>%
  mutate(tOPTV_Zb=as.numeric(pmap(list(CasePTV_Zb, ControlPTV_Zb, 11257, 18501),fisher.or)),
         tPPTV_Zb=as.numeric(pmap(list(CasePTV_Zb, ControlPTV_Zb, 11257, 18501),fisher.pval)))

Jf<-Jf%>%
  replace_na(list(CasePTV_Jf = 0,
                  ControlPTV_Jf = 0)) %>%
  mutate(OPTV_Jf=as.numeric(pmap(list(CasePTV_Jf, ControlPTV_Jf, 11257, 18501),fisher.or)))

###

ggplot(data=tZy,aes(x=as.numeric(tOPTV_Zy),y=as.numeric(OPTV_Zy)))+
  geom_jitter(aes(alpha=0.1))+
  geom_abline()

ggplot(data=tZy,aes(x=-log10(as.numeric(tPPTV_Zy)),y=-log10(as.numeric(PPTV_Zy))))+
  geom_jitter(aes(alpha=0.1))+
  geom_abline()

ggplot(data=tZy,aes(x=-log10(as.numeric(tPPTV_Zy)),y=-log10(as.numeric(PPTV_Zy))))+
  geom_jitter(aes(alpha=0.01,size=100,stroke=20),color="white")+
  geom_abline()


##相关性测试

Zy<-tZy

tZy<-filter(tZy,(!is.na(OPTV_Zy))&(!is.na(tOPTV_Zy)))
tZy<-filter(tZy,(!is.infinite(OPTV_Zy))&(!is.infinite(tOPTV_Zy)))

ggplot(tZy,aes(x=tOPTV_Zy,y=OPTV_Zy))+geom_point()##仅仅是为了验证

cor(tZy$OPTV_Zy,tZy$tOPTV_Zy)
cor(tZy$PPTV_Zy,tZy$tPPTV_Zy)

###最终输出备份

write_csv(Jf,"~/Desktop/INFsaver/fJf.csv")
write_csv(Zb,"~/Desktop/INFsaver/fZb.csv")
write_csv(Zy,"~/Desktop/INFsaver/fZy.csv")

##标准导入格式

Jf<-read_csv("~/Desktop/INFsaver/fJf.csv")
Zb<-read_csv("~/Desktop/INFsaver/fZb.csv")
Zy<-read_csv("~/Desktop/INFsaver/fZy.csv")

##合流JOIN

Fn<-full_join(Jf,Zb,by="Gene")%>%
  full_join(Zy,by="Gene")

##合流备份#1 ##--！谨慎运行！！！！已经注释但不意味着没用

# write_csv(Fn,"~/Desktop/INFsaver/Fn.csv")

#合流备份导入

Fn<-read_csv("~/Desktop/INFsaver/Fn.csv")

#Gene Description join function define
# name_join<-function(DesA,DesB,DesC){
#   # DesA<-as.character(DesA)
#   # DesB<-as.character(DesB)
#   # DesC<-as.character(DesC)
#   if(!is.na(DesA)){
#     DesU<-DesA 
#     return(as.character(DesU))
#   } 
#   if(!is.na(DesB)){
#     DesU<-DesB
#     return(as.character(DesU))
#   } 
#   if(!is.na(DesC)){
#     DesU<-DesC
#     return(as.character(DesU))
#   } 
# }

# #Gene Description Join Test Prepare
# 
# tJf<-fJf
# tZb<-fZb
# tZy<-fZy
# tFn<-full_join(tJf,tZb,by="Gene")%>%
#   full_join(tZy,by="Gene")
# 
# #GENEDESCRIPTIONJOINTEST I
# 
# mutate(tFn,Description_Un=name_join(Description_Jf,Description_Zb,Description_Zy))

#GENE DESCRIPTION UN-lize
# 
# nFn<-mutate(Fn,Description_UN=name_join(Description_Jf,Description_Zb,Description_Zy))%>%
##$有用##select(Gene,CasePTV_Jf,ControlPTV_Jf,CaseMissense_B3_Jf,ControlMissense_B3_Jf,CaseMissense_B23_Jf,ControlMissense_B23_Jf,DenovoPTV_Jf,DenovoMissense_B3_Jf,DenovoMissense_B23_Jf,Pm_Jf,Qm_Jf,OB23_Jf,ORB3_Jf,OMissense_B23_Jf,PMissense_B23_Jf,OMissense_B3_Jf,PMissense_B3_Jf,CaseDenovoPTV_Zb,ControlDenovoPTV_Zb,CaseDenovoMisA_Zb,ControlDenovoMisA_Zb,CaseDenovoMisB_Zb,ControlDenovoMisB_Zb,CasePTVDEN_Zb,ControlPTVDEN_Zb,CasePTVSWE_Zb,ControlPTVSWE_Zb,Transmitted_Zb,Untransmitted_Zb,Qm_Zb,CasePTV_Zb,ControlPTV_Zb,Cases_Zy,Controls_Zy,CasePTV_Zy,ControlPTV_Zy,PPTV_Zy,OPTV_Zy,CaseMissense_Zy,ControlMissense_Zy,PMissense_Zy,OMissense_Zy,tOPTV_Zy,tPPTV_Zy,Description_UN)
# #测试性质
# nFn<-mutate(Fn,Description_UN=name_join(Description_Jf,Description_Zb,Description_Zy))
# Compare<-select(nFn,Description_Jf,Description_Zb,Description_Zy,Description_UN)

#描述合流
Fn<-mutate(Fn,Description_UN=if_else(is.na(Description_Jf),if_else(is.na(Description_Zb),Description_Zy,Description_Zb),Description_Jf))%>%
select(Gene,CasePTV_Jf,ControlPTV_Jf,CaseMissense_B3_Jf,ControlMissense_B3_Jf,CaseMissense_B23_Jf,ControlMissense_B23_Jf,DenovoPTV_Jf,DenovoMissense_B3_Jf,DenovoMissense_B23_Jf,Pm_Jf,Qm_Jf,OB23_Jf,ORB3_Jf,OMissense_B23_Jf,PMissense_B23_Jf,OMissense_B3_Jf,PMissense_B3_Jf,OPTV_Jf,CaseDenovoPTV_Zb,ControlDenovoPTV_Zb,CaseDenovoMisA_Zb,ControlDenovoMisA_Zb,CaseDenovoMisB_Zb,ControlDenovoMisB_Zb,CasePTVDEN_Zb,ControlPTVDEN_Zb,CasePTVSWE_Zb,ControlPTVSWE_Zb,Transmitted_Zb,Untransmitted_Zb,Qm_Zb,CasePTV_Zb,ControlPTV_Zb,OPTV_Zb,PPTV_Zb,Cases_Zy,Controls_Zy,CasePTV_Zy,ControlPTV_Zy,PPTV_Zy,OPTV_Zy,CaseMissense_Zy,ControlMissense_Zy,PMissense_Zy,OMissense_Zy,tOPTV_Zy,tPPTV_Zy,Description_UN)

##########
#筛选数据#
##########

Da<-select(Fn,Gene,Description_UN,OPTV_Jf,Pm_Jf,OPTV_Zb,PPTV_Zb,OPTV_Zy,PPTV_Zy)##半成品没有用管道
Da<-mutate(Da,PPTV_Jf=Pm_Jf)
Ad<-select(Da,Gene,Description_UN,OPTV_Jf,PPTV_Jf,OPTV_Zb,PPTV_Zb,OPTV_Zy,PPTV_Zy)

###最终输出备份

write_csv(Ad,"~/Desktop/INFsaver/Ad.csv")

##最终导入

Ad<-read_csv("~/Desktop/INFsaver/Ad.csv")

library(qqman)

# type  Case  Control   SUM
# JF    24248 97322   121570
# ZB    11257 18501   29758
# ZY    13922 14422   28344

p<-select(Ad,Gene,PPTV_Jf,PPTV_Zb,PPTV_Zy)%>%
  mutate(wJf=348.6689,wZb=160.493,wZy=168.3568)

sp<-p%>%
  filter((!(is.na(PPTV_Jf)))&(!(is.na(PPTV_Zb)))&(!(is.na(PPTV_Zy))))%>%
  rowwise()%>%
  mutate(pvalues=list(c(PPTV_Jf,PPTV_Zb,PPTV_Zy)),weights=list(c(wJf,wZb,wZy)))%>%
  ungroup()%>%
  mutate(sumz_o=map2(pvalues,weights,sumz))%>%
  rowwise()%>%
  mutate(PPTV_UN=sumz_o$p)%>%
  select(-sumz_o,-wJf,-wZb,-wZy)

sp<-select(sp,-pvalues,-weights)
sp<-mutate(sp,PPTV_UN=PPTV_UN[,1])

#sp备份

write_csv(sp,"~/Desktop/INFsaver/Sp.csv")
GS<-read_csv("~/Desktop/INFsaver/Gs.csv")


#####
#选择Gene,pLI,chr

GS<-mutate(GS,chrs=as.numeric(str_replace(chr,"chr","")))%>%
  select(ensembl_gene_id,chrs,pLI)%>%
  mutate(chr=chrs)%>%
  select(-chrs)

GS<-mutate(GS,Gene=ensembl_gene_id)%>%
  select(-ensembl_gene_id)

GS<-select(GS,Gene,chr,pLI)

####
#Join_CH
nsp<-sp
tgs<-GS
fnsp<-left_join(nsp,tgs,by="Gene")


##
write_csv(sp,"~/Desktop/INFsaver/Sp.csv")
write_csv(GS2,"~/Desktop/INFsaver/Gs.csv")

##选择最小

sPPTVgenerater<-function(x,y,z,Q){
  if((x<y)&&(x<z)) temp<-x
  if((y<x)&&(y<z)) temp<-y
  if((z<x)&&(z<y)) temp<-z
  
  if(!is.na(Q)) temp<-NA
  if(x==1&&y==1&&z==1) temp<-NA
  
  if(is.numeric(temp)){
    return(temp)
  }else{
    return(NA)
  }
  
}

nsp<-mutate(sp,sPPTV_UN=sPPTVgenerater(PPTV_Jf,PPTV_Zb,PPTV_Zy,PPTV_UN))


