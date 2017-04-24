TimeSeries<-function
(dates, dateformat, data = NULL, tz = "GMT") 
{
  dimensions<-dim(data)
  clss<-lapply(data,class)
  if (!is.null(data)){
    if (!is.null(dimensions[1])){
      if (length(dates)!=length(data[,1]))
        stop("Data and Date lengths differ")
    }else{
      if (length(dates)!=length(data))
        stop("Data and Date lengths differ")
    }
  }	
  dates <- (strptime(paste(dates), dateformat, tz = tz))
  minute <- minute(dates)
  hour <- hour(dates)
  day <- day(dates)
  week <- week(dates)
  month <- month(dates)
  year <- year(dates)
  if (is.null(data)) {
    results <- data.frame(dates, minute, hour, day, week, 
                          month, year)
  }
  else {
    results <- data.frame(dates, minute, hour, day, week, 
                          month, year, data)
  }
  return(results)
}

daysAgg<-function(data,process,multiple=NULL,na.rm=FALSE){
		if(is.null(multiple)){
		multiple=1
	}
		if(multiple==1){
		day<-aggregate(data[,8:length(data)],list(day=data$day,month=data$month,year=data$year),process,na.rm=na.rm)
	days<- ymd(paste(day$year, day$month, day$day))
	data2<-data.frame(date=days,data=day[,4:length(day)])
	names(data2)<-c("Date",names(data[8:length(data)]))
	return(data2)
		}
	temp<-data
	day<-aggregate(list(data[,8:length(data)],count=1),list(day=data$day,month=data$month,year=data$year),process,na.rm=na.rm)
	days<- ymd(paste(day$year, day$month, day$day))
	data<-data.frame(date=days,day[,5:length(day)-1],count=day[length(day)])
	days=paste(multiple,"days")
	all.dates<-seq.Date(as.Date(data$date[1]),as.Date(data$date[length(data[,1])]),by="day")
	dates<-data.frame(date=all.dates)
	aggreGated<-merge(dates,data,by="date",all.x=TRUE)
	aggreGated$date<-rep(seq.Date(as.Date(data$date[1]),as.Date(data$date[length(data[,1])]),by=days),each=multiple,length=length(all.dates))
#	data<-subset(aggreGated,!is.na(count)) 
	results<-aggregate(list(aggreGated[2:length(aggreGated)]),list(date=aggreGated$date),process,na.rm=TRUE)
	results<-subset(results,results$count!=0)
	results<-results[,-length(results)]
	names(results)<-c("Date",names(temp[8:length(temp)]))
	return(results)
}


hoursAgg<-function 
(data, process, multiple = 1, na.rm = FALSE, tz = "GMT") 
{
	gap <- as.numeric(difftime(data$dates,data$dates[1],tz="GMT", units = "hours"))
	agg.gap<-gap-(gap%%multiple)
	data$dates<-data$dates[1] + agg.gap * 60 * 60
	data <- TimeSeries(data$dates, "%Y-%m-%d %H:%M:%S",
                     data[, 8:length(data)], tz = tz)
	result <- aggregate(data[, 8:length(data)], 
                      list(day = data$day,month = data$month, 
                           year = data$year, hour = data$hour,
                           minute=data$minute),
                      process, na.rm = na.rm)
  data2 <- data.frame(Date = strptime(
    paste(result$minute, result$hour,result$day, result$month, result$year),
    "%M %H %d %m %Y",tz = tz), result[, 6:length(result)])
  sorted <- data2[order(data2$Date), ]
  final <- data.frame(Date = as.POSIXlt(sorted$Date), 
                      sorted[, 2:length(sorted)])
  names(final) <- c("Date", names(data[8:length(data)]))
  return(final)
}

monthsAgg<-function(data,
		    process,
	    	    multiple=NULL,
		    na.rm=FALSE){
	
	if(is.null(multiple)){
		multiple=1
	}	
	
	d.cols<-length(data)
	month<-aggregate(data[,8:d.cols],list(month=data$month,year=data$year),process,na.rm=na.rm)
	data.cols<-length(month)

	
	if(multiple>1){

month.gap<-month[,1]
for(i in 1:length(month[,1])){
	month.gap[i]=(month[i,2]%%month[1,2])*12+month[i,1]}
month.gap<-month.gap-month.gap%%multiple
month.gap<-month.gap%%12+1

year<-month[,2]
if(sum(month.gap)==length(month.gap)){
year<-year[1]+(year-year[1])-(year-year[1])%%(multiple/12)
}else{

for(i in 2:length(month.gap)){
	if(month.gap[i]==month.gap[i-1])
		year[i]=year[i-1]
	else
		year[i]=month[i,2]
		}
	}
		date=strptime(paste(01,month.gap,year),"%d %m %Y")
		results<-data.frame(date,data=month[,3:data.cols])
		final<-aggregate(results[,2:length(results)],list(date=results$date),process,na.rm=na.rm)
		names(final)<-c("Date",names(data[8:length(data)]))
		return(final)
}
	else{
		date=strptime(paste(01,month$month,month$year),"%d %m %Y")
		results<-data.frame(date,data=month[,3:data.cols])
		names(results)<-c("Date",names(data[8:length(data)]))
		return(results)
	}
}

yearsAgg<-function(data,
		    process,
	    	    multiple=NULL,
		    na.rm=FALSE,
		    from.first.obs=TRUE){
	
	if(is.null(multiple)){
		multiple=1
	}
	#Test the last gap to make sure it's complete.
	start<-min(data$year)
	end<-max(data$year)
	if((end-start)%%multiple!=0)warning("Last Gap does not contain ",multiple," years. There is only ",((end-start)%%multiple)," year(s) in this gap.",call.=FALSE)
	
	if(from.first.obs==TRUE){
	years<-(as.numeric(difftime(data$date,data$date[1],units=c("hours")))-as.numeric(difftime(data$date,data$date[1],units=c("hours")))%%8765.81277)/8765.81277
	#8765.81277 hours in each year.
	years.again<-years-years%%multiple
	data$year<-data$year[1]+years.again
	date.1<-as.Date(strptime(paste(data$day[1],data$month[1],data$year),"%d %m %Y"))
	new.1<-data.frame(date=date.1,data[,8:length(data)])
	result<-aggregate(new.1[,2:length(new.1)],list(date=new.1$date),process,na.rm=na.rm)
	names(result)<-c("Date",names(data[8:length(data)]))
	}else{
		years<-data$year-data$year[1]
		years.again<-years-years%%multiple
	data$year<-data$year[1]+years.again
	date.1<-as.Date(strptime(paste(1,1,data$year),"%d %m %Y"))
	new.1<-data.frame(date=date.1,data[,8:length(data)])
	result<-aggregate(new.1[,2:length(new.1)],list(date=new.1$date),process,na.rm=na.rm)
	names(result)<-c("Date",names(data[8:length(data)]))
		}
	

	
		#Everything that need to be done to calculate the annual agg.
	return(result)		
		}


TimeSeries2zoo <- function(x, ...) {
	stopifnot(require(zoo))
	zoo(x[-seq(7)], x$dates)
}

zoo2TimeSeries <- function(x, ...) {
	stopifnot(require(zoo))
	stopifnot(inherits(time(x), "Date") || inherits(time(x), "POSIXt"))
	fmt <- if (inherits(time(x), "Date")) "%Y-%m-%d" else "%Y-%m-%d %H:%M:%S"
	if (length(dim(x)) < 2) {
		DF <- data.frame(coredata(x))
		names(DF) <- deparse(substitute(x))
		TimeSeries(time(x), fmt, DF)
	} else TimeSeries(time(x), fmt, as.data.frame(coredata(x)))
}

