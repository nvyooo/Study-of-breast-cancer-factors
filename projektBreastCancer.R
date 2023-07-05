install.packages("bnlearn")
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(c("graph", "RBGL", "Rgraphviz"))
install.packages("gRain")
install.packages("pcalg")
library(bnlearn)
library(Rgraphviz)
library(gRain)
library(pcalg)

#-----Wczytanie danych-----
dane<- read.csv2("C:/Users/theal/Desktop/ProjektBreastCancer/Breast_Cancer.csv", sep = ",")
View(dane)

#-----zmiana na factor-----
str(dane)
dane <- lapply(dane,factor)
dane <- as.data.frame(dane)
str(dane)

#-----algorytm hc-----
siec_hc<-hc(dane)

#graf
graphviz.plot(siec_hc)

# struktura sieci
modelstring(siec_hc)
nodes(siec_hc)


#Badanie siły łuków sieci
arc.strength(siec_hc, data=dane)

#dopasowanie struktury do danych
score(siec_hc, data=dane, type="bic")

model_hc<-bn.fit(siec_hc, dane)

graphviz.chart(model_hc, type = "barprob", grid = TRUE,
               bar.col = "darkgreen", strip.bg = "lightskyblue")


#-----algorytm tabu-----
siec_tabu <- tabu(dane)
print(siec_tabu)

# graf
graphviz.plot(siec_tabu)

# badanie siły łuków sieci
arc.strength(siec_tabu, data = dane)

# struktura sieci
modelstring(siec_tabu)

# dopasowanie struktury do danych
score(siec_tabu, data=dane, type="bic")

#dopasowanie parametrów sieci
model_tabu <- bn.fit(siec_tabu, dane)

#-----algorytm aracne I -----

siec_aracne <- aracne(dane)
print(siec_aracne)

graphviz.plot(siec_aracne)


# testy chi kwadrat na sprawdzenie zależności w grafie
ci.test("Age","Race", data=dane, test="x2")#brak podstaw
ci.test("Age","differentiate", data=dane, test="x2")#brak podstaw
ci.test("Age","Grade", data=dane, test="x2") #brak podstaw
ci.test("Race", "differentiate", data=dane, test="x2") #brak podstaw
ci.test("Race", "Grade", data=dane, test="x2") #brak podstaw
ci.test("Grade", "differentiate", data=dane, test="x2") #brak podstaw
ci.test("Marital.Status", "Age", data=dane, test="x2") #brak podstaw
ci.test("Age", "T.Stage", data=dane, test="x2") #odrzucamy
ci.test("Age", "N.Stage", data=dane, test="x2") #odrzucamy
ci.test("Age", "X6th.Stage", data=dane, test="x2")#odrzucamy
ci.test("Race","Tumor.Size", data=dane, test="x2") #odrzucamy
ci.test("Age","Tumor.Size", data=dane, test="x2")#odrzucamy
ci.test("T.Stage","N.Stage", data=dane, test="x2") #brak podstaw


# Definiowanie łuków w grafie nieskierowanym
siec_aracne <- set.arc(siec_aracne, "Age","Race")
siec_aracne <- set.arc(siec_aracne, "Age","differentiate")
siec_aracne <- set.arc(siec_aracne, "Age","Grade")
siec_aracne <- set.arc(siec_aracne, "Race", "differentiate")
siec_aracne <- set.arc(siec_aracne, "Race", "Grade")
siec_aracne <- set.arc(siec_aracne, "Marital.Status","N.Stage")
siec_aracne <- set.arc(siec_aracne, "Marital.Status","Tumor.Size")
siec_aracne <- set.arc(siec_aracne, "Estrogen.Status","Tumor.Size")
siec_aracne <- set.arc(siec_aracne, "A.Stage","Tumor.Size")

# Tworzenie pustego grafu skierowanego
siec_aracne_directed <- empty.graph(nodes = nodes(siec_aracne))

# Dodawanie łuków na podstawie grafu nieskierowanego
arcs <- arcs(siec_aracne)
for (i in 1:nrow(arcs)) {
  from <- arcs[i, "from"]
  to <- arcs[i, "to"]
  siec_aracne_directed <- set.arc(siec_aracne_directed, from, to)
}

# Wyświetlanie węzłów w grafie skierowanym
nodes <- nodes(siec_aracne_directed)
print(nodes)

# Wyświetlanie łuków w grafie skierowanym
arcs <- arcs(siec_aracne_directed)
print(arcs)

# Struktura sieci
modelstring(siec_aracne_directed)

# graf
graphviz.plot(siec_aracne_directed)

# dopasowanie struktury do danych
score(siec_aracne_directed,data=dane,type="bic")

#dopasowanie parametrów sieci
model_aracne <- bn.fit(siec_aracne_directed, dane)

#dopasowanie parametrów sieci
model_tabu <- bn.fit(siec_tabu, dane)

#------Algorytm GS - Grow-Shrink------ 
# -------------------------
siec_gs <- gs(dane)
# print(siec_gs)

# graf
graphviz.plot(siec_gs)  


# Definiowanie łuków w grafie nieskierowanym
siec_gs <- set.arc(siec_gs, "Age","Tumor.Size")
siec_gs <- set.arc(siec_gs, "Age", "T.Stage")
siec_gs <- set.arc(siec_gs, "Race","Tumor.Size")

# Tworzenie pustego grafu skierowanego
siec_gs_directed <- empty.graph(nodes = nodes(siec_gs))

# Dodawanie łuków na podstawie grafu nieskierowanego
gs <- arcs(siec_gs)
for (i in 1:nrow(arcs)) {
  from <- arcs[i, "from"]
  to <- arcs[i, "to"]
  siec_gs_directed <- set.arc(siec_gs_directed, from, to)
}

# dopasowanie struktury do danych
score(siec_gs_directed,data=dane,type="bic") 


# struktura sieci
modelstring(siec_gs_directed)

# graf
graphviz.plot(siec_gs_directed)     

#dopasowanie parametrów sieci
model_gs <- bn.fit(siec_gs_directed, dane)


# -------------------------------------
# pc.stable i ręczne definiowanie łuków
# -------------------------------------
siec <- pc.stable(dane)

# graf
graphviz.plot(siec)

# Tworzenie grafu nieskierowanego
siec <- empty.graph(nodes = colnames(dane))
siec <- set.arc(siec, "Age", "Race")
siec <- set.arc(siec, "Race", "Grade")
siec <- set.arc(siec, "Estrogen.Status", "Tumor.Size")


# badanie siły łuków sieci
arc.strength(siec, data = dane)

# dopasowanie struktury do danych
score(siec,data=dane,type="bic") 



# ---------------------
# whitelist i blacklist
# ---------------------
sieclista <- hc(dane, whitelist=matrix(c(c("Race", "Grade"),
                                         c("Marital.Status","Tumor.Size")),
                                       ncol = 2, byrow = TRUE),
                blacklist = matrix(c(c("Age", "Race"),
                                     c("Race", "differentiate")),
                                   ncol = 2, byrow = TRUE))

# graf
graphviz.plot(sieclista) 

# dopasowanie struktury do danych
score(sieclista,data=dane,type="bic")



# MLE - maximum likelihood estimation
bn <- empty.graph(colnames(dane))
bn.fit <- bn.fit(bn, dane, method="mle")
print(bn.fit)


#model_hc <- bn.fit(siec_hc,dane)

#---------prawdopodobienstwa-------

# Przekształcenie zmiennych Factor na numeryczne
dane$Age <- as.numeric(dane$Age)
dane$Tumor.Size <- as.numeric(dane$Tumor.Size)
dane$Estrogen.Status <- as.numeric(dane$Estrogen.Status)
dane$Race <- as.numeric(dane$Race)
dane$Grade <- as.numeric(dane$Grade)
dane$Status <- as.numeric(dane$Status)

str(dane)

# Prawdopodobieństwo, że losowo wybrany pacjent miał wysoki Tumor Size
model_hc$Tumor.Size
junction <-compile(as.grain(model_hc))
querygrain(junction, nodes = "Estrogen.Status")$Estrogen.Status

querygrain(junction, nodes = "Tumor.Size")$Tumor.Size

# Prawdopodobieństwo, że losowo wybrany pacjent, który ma powyżej 30 lat,
# będzie miał wysoki Tumor Size, niski Estrogen Status oraz będzie biały
numerator <- sum(dane$Age > 30 & dane$Tumor.Size > 25 & dane$Estrogen.Status == 2 & dane$Race == 1)
denominator <- sum(dane$Age > 30 & dane$Tumor.Size > 0 & dane$Estrogen.Status == 2 & dane$Race == 1)
result <- numerator / denominator
result

# Prawdopodobieństwo, że losowo wybrany pacjent miał Grade 3 i wysoki Estrogen Status,
# jeśli umarł
sum(dane$Grade == 3 & dane$Estrogen.Status == 1 & dane$Status == 2) / sum(dane$Status == 1)
model_hc$Grade


