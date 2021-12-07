#' Cross-national data on remittances and protest
#'
#' A data set to replicate the findings of Escrib\`{a}-Folch,
#' Meseguer, and Wright (2018). Data and data descriptions are
#' from that paper's replication data, available at  
#' \doi{10.7910/DVN/TVZQG6}
#'
#'
#' @docType data
#' @name remit
#' @format  A data frame with 2429 observations and 14 variables:
#' \describe{
#'  \item{Protest}{standardized measure of latent protest from Chenoweth et al. (2014)}
#'  \item{remit}{natural log of the 2-year lagged moving average of total remittances received in constant US dollars}
#'  \item{dict}{binary indicator of autocracy or democracy from	Geddes, Wright, and Frantz (2014)}
#'  \item{l1gdp}{natural log of one-period lagged gdp per capita}
#'  \item{l1pop}{natural log of one-period lag of population}
#'  \item{l1nbr5}{lagged mean latent level of protest in countries with capital cities within 4000km of the target country's capital}
#'  \item{l12gr}{two-year lagged moving average of GDP per capita growth (in percent)}
#'  \item{l1migr}{natural log of lagged net migration in millions}
#'  \item{elec3}{indicator for multiparty election in that year, year prior, or year after}
#'  \item{cowcode}{country code from correlates of war dataset}
#'  \item{period}{six ordinal time periods}
#'  \item{caseid}{numerical code for autocratic regime case name}
#'  \item{year}{year}
#' }
#' @source \doi{10.7910/DVN/TVZQG6}
#' @references Escrib\`{a}-Folch, A., Meseguer, C. and Wright, J. (2018), Remittances and Protest in Dictatorships. American Journal of Political Science, 62: 889-904. \doi{10.1111/ajps.12382}
#'
#' Wright, Joseph, 2018, "Replication Data for: Remittances and Protest in Dictatorships", \doi{10.7910/DVN/TVZQG6}, Harvard Dataverse, V1, UNF:6:IE6OqUb3EB5AIDYKI28mgA== [fileUNF]
#' 
"remit"


#' Data on the direct primary in US congressional elections
#'
#' A data set on the presence of the direct primary in U.S.
#' congressional elections and the vote shares for the Democratic,
#' Republican, and third parties. Based on ICPSR Study 6985
#'
#' @docType data
#' @name primary
#' @format  A data frame with 1164 observations and the following 7
#'   variables:
#' \describe{
#'   \item{state}{name of the state}
#'   \item{year}{year of the congressional election}
#'   \item{dem_share}{percentage of the total vote cast for the
#'   Democratic candidate, 0-100}
#'   \item{rep_share}{percentage of the total vote cast for the
#'   Republican candidate, 0-100}
#'   \item{other_share}{percentage of the total vote cast for other
#'   parties, 0-100}
#'   \item{primary}{binary variable indicating if the state had the
#'   direct primary (=1) or not (=0)}
#'   \item{south}{binary variable indicating if the state is in the
#'   South (=1) or not (=0)}
#' }
#' @source \url{https://www.icpsr.umich.edu/icpsrweb/ICPSR/studies/6895}
#' @references David, Paul T., and Claggett, William. Party Strength
#'   in the United States: 1872-1996. Ann Arbor, MI: Inter-university
#'   Consortium for Political and Social Research [distributor],
#'   2008-09-10.  https://doi.org/10.3886/ICPSR06895.v1
#'
#' 
"primary"
