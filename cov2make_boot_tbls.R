nextclade_path <- 'PATH TO NEXTCLADE OUTPUT IN *.tsv FORMAT AFTER FILTRATION, *.Nfiltered.tsv'
meta_path <- 'PATH TO METADATA IN GISAID *.tsv FORMAT'
write_path = 'PATH TO WRITE gisaid_nc_mrg'
plus1_path = 'PATH TO cov2plus1_ptrns output, "from-to" table'
world_boot_tbl_path = 'PATH TO world boot table'
rus_boot_tbl_path = 'PATH TO Russian boot table'
rus_boot_tbl_with_fd_path  = 'PATH TO Russian boot table with FD'

nc_gisaid <- read.table(nextclade_path, sep = '\t', header = TRUE, quote = "", comment.char = "")
gisaid_metadata <- read.table(meta_path, sep = '\t', header = TRUE, quote = "", comment.char = "")


nc_gisaid$substitutions[nc_gisaid$substitutions == ''] <- 'root'

if(is.null(nc_gisaid$seqName_short)) nc_gisaid$seqName_short <- gsub('\\|.*', '', nc_gisaid$seqName)

gisaid_nc_mrg <- merge(nc_gisaid, gisaid_metadata, by.x = 'seqName_short', by.y = 'Virus.name', all.x = TRUE)
gisaid_nc_mrg$Location <- gsub(' ', '', gisaid_nc_mrg$Location)
gisaid_nc_mrg$Location <- gsub(' |-|–', '', gisaid_nc_mrg$Location)
gisaid_nc_mrg$country <- sapply(gisaid_nc_mrg$Location, function(x) strsplit(x, "/")[[1]][2])
gisaid_nc_mrg$country <- sapply(gisaid_nc_mrg$country, function(x) gsub("'", "", x))
gisaid_nc_mrg$country <- sapply(gisaid_nc_mrg$country, function(x) gsub("\\.", "", x))



regs <- read.table('countriy_regions.tsv', sep = '\t', header = TRUE, quote = '')



gisaid_nc_mrg$subregion <- sapply(gisaid_nc_mrg$country, function(x) ifelse(x %in% regs$Country.or.Area, regs$Geographical.subregion[regs$Country.or.Area == x], NA))


# SEPARATE CONTRIES WITH BIG NUMBER OF GENOMES
# big_cntrs <- c('India(SouthernAsia)',
#                  'Spain(SouthernEurope)',
#                  'Italy(SouthernEurope)',
#                  'Japan(EasternAsia)',
#                  'SouthKorea(EasternAsia)',
#                  'Brazil(SouthAmerica)',
#                  'Turkey(WesternAsia)',
#                  'Israel(WesternAsia)',
#                  'Germany(WesternEurope)',
#                  'France(WesternEurope)',
#                  'Poland(EasternEurope)',
#                  'Russia(EasternEurope)',
#                  'USA(NorthernAmerica)',
#                  'UnitedKingdom(NorthernEurope)',
#                  'Denmark(NorthernEurope)',
#                  'Mexico(CentralAmerica)',
#                  'Singapore(SouthEasternAsia)')
# gisaid_nc_mrg$subregion <- sapply(gisaid_nc_mrg$country, function(x){
#   if(x %in% big_cntrs){
#     return(big_cntrs[which(big_cntrs == x)])
#   } else {
#     ifelse(x %in% regs$Country.or.Area, regs$Geographical.subregion[regs$Country.or.Area == x], NA)
#   }
# })

# SEPARATE CONTRIES WITH BIG NUMBER OF GENOMES AND COUNTRIES FROM SOUTH-EASTERN ASIA
# big_cntrs.SEA <- c('India(SouthernAsia)',
#                     'Spain(SouthernEurope)',
#                     'Italy(SouthernEurope)',
#                     'Japan(EasternAsia)',
#                     'SouthKorea(EasternAsia)',
#                     'Brazil(SouthAmerica)',
#                     'Turkey(WesternAsia)',
#                     'Israel(WesternAsia)',
#                     'Germany(WesternEurope)',
#                     'France(WesternEurope)',
#                     'Poland(EasternEurope)',
#                     'Russia(EasternEurope)',
#                     'USA(NorthernAmerica)',
#                     'UnitedKingdom(NorthernEurope)',
#                     'Denmark(NorthernEurope)',
#                     'Mexico(CentralAmerica)',
#                     'Singapore(SouthEasternAsia)',
#                     "Indonesia(SouthEasternAsia)",
#                     "Malaysia(SouthEasternAsia)",
#                     "Thailand(SouthEasternAsia)",
#                     "Philippines(SouthEasternAsia)",
#                     "Brunei(SouthEasternAsia)",
#                     "Vietnam(SouthEasternAsia)")
# gisaid_nc_mrg$subregion <- sapply(gisaid_nc_mrg$country, function(x){
#   if(x %in% big_cntrs.SEA){
#     return(big_cntrs.SEA[which(big_cntrs.SEA == x)])
#   } else {
#     ifelse(x %in% regs$Country.or.Area, regs$Geographical.subregion[regs$Country.or.Area == x], NA)
#   }
# })

# SAVE INTERMEDIATE RESULTS
write.table(gisaid_nc_mrg, write_path, sep = '\t', row.names = FALSE, quote = FALSE)


p1m1_df <- read.table(plus1_path, sep = '\t', header = TRUE, quote = '')



p1m1_ptrns <- c(p1m1_df$from, unss(p1m1_df$to, split = ';')) %>% unique
dup_ptrns <- gisaid_nc_mrg$substitutions[duplicated(gisaid_nc_mrg$substitutions)] %>% unique
world_ptrns <- c(p1m1_ptrns, dup_ptrns)
world_boot_tbl <- gisaid_nc_mrg[, c('substitutions', 'subregion')]


setDT(world_boot_tbl)[!substitutions %in% world_ptrns, substitutions := NA]


write.table(world_boot_tbl, world_boot_tbl_path, sep = '\t', row.names = FALSE, quote = FALSE)




############################################### BOOT TABLE FOR RUSSIA
nc_gisaid_rus <-  nc_gisaid[grepl('Russia', nc_gisaid$seqName), ]

gisaid_nc_rus_mrg <- gisaid_nc_mrg[gisaid_nc_mrg$country == 'Russia' & !is.na(gisaid_nc_mrg$country), ]
gisaid_nc_rus_mrg$region <- sapply(gisaid_nc_rus_mrg$Location, function(x) gsub('/.*', '', gsub('Europe/Russia/', '', x)))


######################################RUS_P1M1
p1m1_ptrns_rus <- p1m1_ptrns[p1m1_ptrns %in% gisaid_nc_rus_mrg$substitutions] %>% unique

p1m1_df_rus <- p1m1_df[p1m1_df$from %in% p1m1_ptrns_rus, ]
p1m1_df_rus <- p1m1_df_rus[sapply(p1m1_df_rus$to, function(x) as.logical(sum(unss(x, split = ';') %in% p1m1_ptrns_rus))), ]
p1m1_df_rus$to <- sapply(p1m1_df_rus$to, function(x){
  splt <- unss(x, split = ';')
  paste0(splt[(splt %in% p1m1_ptrns_rus)], collapse = ';')
})

p1m1_df_rus <- p1m1_df_rus[p1m1_df_rus$to != '', ]
p1m1_ptrns_rus <- c(p1m1_df_rus$from, unss(p1m1_df_rus$to, split = ';')) %>% unique
# write.table(p1m1_df_rus, 'PATH_TO_SAVE', sep = '\t', row.names = FALSE, quote = FALSE)
######################################



library(stringi)

fix_reg <- function(vect){
  ptrn <- c("Moscowregion|MoscowRegion","Moscow","VologdaRegion|Volgogradregion|VolgogradRegion|VolgogradOblast|Volgograd","Kalugaregion|Kaluga","SaintPetersburg|St.Petersburg",
            "YamaloNenetsAutonomousOkrug|YamaloNenets","ArkhangelskRegion|ArkhangelskOblast|Arkhangelsk","PenzaRegion|Penza","Lipetskregion.*|Lipetsk","Amurregion|Amur", 
            "RepublicofKalmykia|Kalmykia","NovosibirskRegion|NovosibirskOblast|Novosibirsk","RyazanRegion|Ryazan","MariElRepublic|MariEl",
            "RepublicofBashkortostan|Bashkortostan","IvanovoRegion|Ivanovo","UlyanovskRegion|UlianovskRegion|Ulyanovsk","Tularegion|Tula",
            "VoronezhRegion|VoronezhOblast|Voronezh","RepublicofMordovia|MordoviaRepublic|MordoviyaRepublic|Mordovia","KamchatkaKrai|KamchatkaTerritory|Kamchatka","Tambovregion|Tambov","PrimorskyKrai|PrimorskiyRegion|PrimorskyKra|Primorsky",
            "Smolenskregion|SmolenskOblast|Smolensk","TulaRegion","Magadanregion","TyumenRegion|Tumenregion|Tyumen","RepublicofKarelia|Kareliarepublic|Karelia",
            "LeningradRegion|LeningradOblast|Leningrad","MagadanRegion|Magadan","SaratovRegion|Saratov",
            "Yaroslavlregion|Yaroslavl", "KaliningradRegion|Kalinigradregion|Kaliningrad", "ZabaykalskyKrai|Zabaykalsky", "KrasnodarKrai|Krasnodar", "Vologda", "Kurskregion|KurskRegion|Kursk", "KemerovoRegion|Kemerovo",
            "NizhnyNovgorodRegion|NizhnyNovgorodOblast|NizhnyNovgorod", "NovgorodRegion|NovgorodOblast|Novgorod", "KhantyMansiAutonomousOkrugYugra|KhantyMansiAutonomousOkrug|KhantyMansiAutonomousArea|KhantyMansi",
            "Tverregion|Tver", "Vladimirregion|VladimirOblast|Vladimir", "RepublicofBuryatia|Buryatia", "Kostromaregion|Kostroma", "Rostovregion|Rostov", "KarachayCherkess", "PermKrai|PermRegion|Perm", "ChukotkaA.O.|ChukotkaAO|ChukotkaAutonomousOkrug|Chukotka", "AstrakhanRegion|Astrakhan",
            "RepublicofTatarstan|Tatarstan", "OrenburgRegion|Orenburg", "MurmanskRegion|Murmansk", "Orlovregion|Oryolregionregion|OryolRegion|Oryol|Orel", "Bryanskregion|Bryansk", "KabardinoBalkarianRepublic|KabardinoBalkariaRepublic|RepublicofKabardinoBalkaria|KabardinoBalkaria",
            "RepublicofDagestan|Dagestan", "PskovRegion|Pskov", "StavropolskiyKrai|StavropolKrai", "KirovRegion|Kirov", "ChelyabinskRegion|Chelyabinsk",
            "BelgorodRegion|Belgorod", "ChechenRepublic", "Tomskregion", "Omskregion|Omsk", "Sevastopol", "RepublicofIngushetia|IngushRepublic|Ingushetia",
            "JewishAutonomousOblast", "UdmurtianRepublic|UdmurtRepublic", "TyvaRepublic", "SverdlovskRegion|Sverdlovsk", "RepublicofNorthOssetiaAlania", "RepublicofSakha\\(Yakutia\\)|Sakha\\(Yakutia\\)Republic|SakhaRegion|Yakutia|Yakutiya|Sakha\\(Yakutia\\)",
            "Samararegion|Samara", "SakhalinOblast|SakhalinRegion|Sakhalin", "KomiRepublic|Komi", "NenetsAutonomousOkrug|NenetsAutonomousArea", "KrasnoyarskKrai|Krasnoyarsk",
            "RepublicofKhakassia|Khakassia", "KhabarovskKrai|Khabarovsk", "KarachayevoCircassianRepublic", "KurganRegion", "IrkutskRegion|Irkutskregion|Irkutsk",
            "ChuvashRepublic|Chuvashia", "AltaiskiyKrai|AltaiKrai|AltaiK|Altairegion|Altai", "RepublicofAdygea", "Crimea")
  rplcmnt <- c("MOS","MOW","VLG","KLU","SPE",
               "YAN","ARK","PNZ","LIP","AMU",
               "KL","NVS","RYA","ME",
               "BA","IVA","ULY","TUL",
               "VOR","MO","KAM","TAM","PRI",
               "SMO","TUL","MAG","TYU","KR",
               "LEN","MAG","SAR",
               "YAR", "KGD", "ZAB", "KDA", "VLG", "KRS", "KEM",
               "NIZ", "NGR", "KHM",
               "TVE", "VLA", "BU", "KOS", "ROS", "KC", "PER", "CHU", "AST",
               "TA", "ORE", "MUR", "ORL", "BRY", "KB",
               "DA", "PSK", "STA", "KIR", "CHE",
               "BEL", "CE", "TOM", "OMS", "40", "IN",
               "YEV", "UD", "TY", "SVE", "SE", "SA",
               "SAM", "SAK", "KO", "NEN", "KYA",
               "KK", "KHA", "KC", "KGN", "IRK",
               "CU", "ALT", "AD", "43")
  stri_replace_all_regex(vect, pattern = ptrn, replacement = rplcmnt, vectorize=FALSE, case_insensitive = TRUE)
}

gisaid_nc_rus_mrg$region <- fix_reg(gisaid_nc_rus_mrg$region)
dup_ptrns_rus <- gisaid_nc_rus_mrg$substitutions[duplicated(gisaid_nc_rus_mrg$substitutions)] %>% unique


rus_boot_tbl <- gisaid_nc_rus_mrg[, c('substitutions', 'region')]
rus_boot_tbl$substitutions <- sapply(rus_boot_tbl$substitutions, function(x) ifelse(x %in% c(p1m1_ptrns_rus.v2, dup_ptrns_rus), x, NA))

write.table(rus_boot_tbl, rus_boot_tbl_path, sep = '\t', row.names = FALSE, quote = FALSE)
############################################### BOOT TABLE FOR RUSSIA - END



##########################################################################################################################WORLD BOOT TBL WITH RUS FEDERAL DISTRICTS

iso_to_federal_district <- function(iso_codes) {

  iso_codes <- as.character(iso_codes)
  iso_codes <- toupper(iso_codes)
  
  patterns <- c(
    "^MOW$",
    # Central Federal District (without MOW)
    "^BEL$", "^BRY$", "^VLA$", "^VOR$", "^IVA$", "^KLU$", "^KOS$",
    "^KRS$", "^LIP$", "^MOS$","^ORL$", "^RYA$", "^SMO$", "^TAM$",
    "^TVE$", "^TUL$", "^YAR$", 
    # Northwestern Federal District
    "^KR$", "^KO$", "^ARK$", "^VLG$", "^KGD$", "^LEN$", "^MUR$",
    "^NGR$", "^PSK$", "^SPE$", "^NEN$", 
    # Southern Federal District
    "^AD$", "^KL$", "^43$", "^KDA$", "^AST$", "^VGG$", "^ROS$",
    "^40$",
    # North Caucasian Federal District
    "^DA$", "^IN$", "^KB$", "^KC$", "^SE$", "^CE$", "^STA$",
    # Volga Federal District
    "^BA$", "^ME$", "^MO$", "^TA$", "^UD$", "^CU$", "^PER$",
    "^KIR$", "^NIZ$", "^ORE$", "^PNZ$", "^SAM$", "^SAR$", "^ULY$",
    # Ural Federal District
    "^KGN$", "^SVE$", "^TYU$", "^CHE$", "^KHM$", "^YAN$",
    # Siberian Federal District
    "^AL$", "^TY$", "^KK$", "^ALT$", "^KYA$", "^IRK$",
    "^KEM$", "^NVS$", "^OMS$", "^TOM$",
    # Far Eastern Federal District
    "^BU$", "^SA$", "^ZAB$", "^KAM$", "^PRI$", "^KHA$", "^AMU$", "^MAG$", "^SAK$", 
    "^YEV$", "^CHU$"
  )
  

  replacements <- c(
    "Moscow",
    rep("CentralFD", 17),
    rep("NorthwesternFD", 11),
    rep("SouthernFD", 8),
    rep("NorthCaucasianFD", 7),
    rep("VolgaFD", 14),
    rep("UralFD", 6),
    rep("SiberianFD", 10),
    rep("FarEasternFD", 11)
  )

  result <- stringi::stri_replace_all_regex(
    iso_codes, patterns, replacements,
    vectorize_all = FALSE
  )
  
  unchanged <- result == iso_codes
  result[unchanged] <- NA_character_
  
  return(result)
}

gisaid_nc_mrg_with_rus_fd <- gisaid_nc_mrg

gisaid_nc_mrg_with_rus_fd$region_rus <- sapply(gisaid_nc_mrg_with_rus_fd$Location, function(x){ 
  if(grepl('Russia', x)){gsub('/.*', '', gsub('Europe/Russia/', '', x))
  } else NA
})

gisaid_nc_mrg_with_rus_fd$region_rus <- NA
gisaid_nc_mrg_with_rus_fd$region_rus[!is.na(gisaid_nc_mrg_with_rus_fd$region_rus)] <- fix_reg(gisaid_nc_mrg_with_rus_fd$region_rus[!is.na(gisaid_nc_mrg_with_rus_fd$region_rus)])

# quick fix
# gisaid_nc_mrg_with_rus_fd <- gisaid_nc_mrg_with_rus_fd[!grepl('hCoV-19/Russia/un', gisaid_nc_mrg_with_rus_fd$seqName_short), ]


gisaid_nc_mrg_with_rus_fd$rus_fd <- NA
gisaid_nc_mrg_with_rus_fd$rus_fd[!is.na(gisaid_nc_mrg_with_rus_fd$region_rus)] <- iso_to_federal_district(gisaid_nc_mrg_with_rus_fd$region_rus[!is.na(gisaid_nc_mrg_with_rus_fd$region_rus)])

gisaid_nc_mrg_with_rus_fd$subregion_with_fd <- apply(gisaid_nc_mrg_with_rus_fd[, c('subregion', 'rus_fd')], 1, function(x){
  if(x[1] != 'Russia(EasternEurope)') return(x[1])
  if(x[1] == 'Russia(EasternEurope)') return(paste0(x[2], '(Russia)'))
})

gisaid_nc_mrg_with_rus_fd$subregion_with_fd %>% table

# saving interim results
# write.table(gisaid_nc_mrg_with_rus_fd, 'PATH', sep = '\t', row.names = FALSE, quote = FALSE)

world_boot_tbl_with_rus_fd <- gisaid_nc_mrg_with_rus_fd[, c('substitutions', 'subregion_with_fd')]

setDT(world_boot_tbl_with_rus_fd)[!substitutions %in% world_ptrns, substitutions := NA]

write.table(world_boot_tbl_with_rus_fd, rus_boot_tbl_with_fd_path, sep = '\t', row.names = FALSE, quote = FALSE)










# #############################BEFORE AND AFTER 2022 - example
# gisaid_nc_mrg$Date.group <- sapply(gisaid_nc_mrg$Collection.date, function(x){
#   if(as.numeric(substr(x, 1, 4)) %in% (2019:2021)){
#     return("Before")
#   } else if(as.numeric(substr(x, 1, 4)) %in% (2022:2025)){
#     return('After')
#   } else return("Error")
# })
# 
# 
# nc_gisaid_after <-  gisaid_nc_mrg[gisaid_nc_mrg$Date.group == 'After', ]
# nc_gisaid_before <-  gisaid_nc_mrg[gisaid_nc_mrg$Date.group == 'Before', ]
# 
# 
# gisaid_nc_after_mrg <- gisaid_nc_mrg[gisaid_nc_mrg$seqName_short %in% nc_gisaid_after$seqName_short, ]
# gisaid_nc_before_mrg <- gisaid_nc_mrg[gisaid_nc_mrg$seqName_short %in% nc_gisaid_before$seqName_short, ]
# 
# 
# p1m1_ptrns <- c(p1m1_df$from, unss(p1m1_df$to, split = ';')) %>% unique
# p1m1_ptrns_after <- p1m1_ptrns[p1m1_ptrns %in% gisaid_nc_after_mrg$substitutions] %>% unique
# p1m1_ptrns_before <- p1m1_ptrns[p1m1_ptrns %in% gisaid_nc_before_mrg$substitutions] %>% unique
# 
# 
# p1m1_df_after <- p1m1_df[p1m1_df$from %in% p1m1_ptrns_after, ]
# # write.table(p1m1_df_after, 'WHERE_TO_WRITE_AFTER', sep = '\t', row.names = FALSE, quote = FALSE)
# 
# 
# 
# p1m1_df_before <- p1m1_df[p1m1_df$from %in% p1m1_ptrns_before, ]
# # write.table(p1m1_df_before, 'WHERE_TO_WRITE_BEFORE', sep = '\t', row.names = FALSE, quote = FALSE)
# 
# 
# ###############################################################################################AFTER
# p1m1_ptrns_after <- c(p1m1_df_after$from, unss(p1m1_df_after$to, split = ';')) %>% unique
# dup_ptrns_after <- gisaid_nc_after_mrg$substitutions[duplicated(gisaid_nc_after_mrg$substitutions)] %>% unique
# world_ptrns_after <- c(p1m1_ptrns_after, dup_ptrns_after)
# 
# world_boot_tbl.AFTER <- gisaid_nc_after_mrg[, c('substitutions', 'subregion.v2')]
# setDT(world_boot_tbl.AFTER)[!substitutions %in% world_ptrns_after, substitutions := NA]
# write.table(world_boot_tbl.AFTER, 'WHERE_TO_WRITE_AFTER_BOOT_TABLE', sep = '\t', row.names = FALSE, quote = FALSE)
# ##############################################################################################BEFORE
# p1m1_ptrns_before <- c(p1m1_df_before$from, unss(p1m1_df_before$to, split = ';')) %>% unique
# dup_ptrns_before <- gisaid_nc_before_mrg$substitutions[duplicated(gisaid_nc_before_mrg$substitutions)] %>% unique
# world_ptrns_before <- c(p1m1_ptrns_before, dup_ptrns_before)
# 
# world_boot_tbl.BEFORE <- gisaid_nc_before_mrg[, c('substitutions', 'subregion.v2')]
# setDT(world_boot_tbl.BEFORE)[!substitutions %in% world_ptrns_before, substitutions := NA]
# write.table(world_boot_tbl.BEFORE, 'WHERE_TO_WRITE_BEFORE_BOOT_TABLE', sep = '\t', row.names = FALSE, quote = FALSE)




