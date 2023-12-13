yun_data <- read_excel("data/raw/PademosTibetTraitData-Dryad2019.10.25.xlsx", 
                       skip = 1)

temp <- yun_data %>% group_by(Plot_No) %>% summarise(total_cover = sum(`Cover_%`))

austraits <- load_austraits(version = "3.0.2", path = "data/austraits")
my.ausplots.data <- try(get_ausplots(bounding_box = c(125, 140, -40, -10)))
traits_temp <- austraits$traits

length(unique(traits_temp$taxon_name))
temp2 <- my.ausplots.data$veg.PI

length(unique(temp2$standardised_scientific_name))

Appendix_S2_C_sites <- read_csv("Potential_datasets/Guerin et al. 2022/Appendix_S2_C_sites.csv")
Appendix_S2_D_veg <- read_csv("Potential_datasets/Guerin et al. 2022/Appendix_S2_D_veg.csv")
Appendix_S2_A_traits <- read_csv("Potential_datasets/Guerin et al. 2022/Appendix_S2_A_traits.csv")

relative_abundances_genuslvl <- calculate_abundance(veg_genuslvl, traits_genuslvl)
relative_abundances_familylvl <- calculate_abundance(veg_familylvl, traits_familylvl)