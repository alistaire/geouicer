#library(raster)
#library(rgdal)
#library(sf)
#library(sp)
#library(dplyr)
#library(geojsonsf)
#library(jsonlite)
#library(aws.s3)

RADIUS_OF_EARTH  = 6378137

#' calculates the radius of the earth given a latitude
#'
#'  Calculates the radius at a given latitude assuming a
#'  sphere with the radius of the earth
#'
#'  @param latitude latitude in degrees
#'  @return radius of the earth at latitude in meters
#'
#'  @examples
#'  radius_at_latitude(-56)
#'
radius_at_latitude <- function(latitude){
  lrad = latitude/360 * 2*pi
  result = cos(lrad)*RADIUS_OF_EARTH
  return (result)
}


#' Calculates steps for incremental buffering
#'
#' Calculates a set of steps to to buffer up to a given distance.
#' Can be more efficient than buffering in one go if there are a lot
#' of features to buffer and so a lot interaction between features
#' likely to be created when buffers are unioned. Heuristic function
#' likely to be better approaches. Rule of Thumb to increment in 1s and 3
#' is another approach e.g. 1, 3, 10, 30, 100, 300, 1000, 3000, ....
#' might also work
#'
#' @param mn - smallest step to start from
#' @param mx - overall buffer distance required
#' @param num_step - number of increments
#' @return matrix of steps that sum mx
#'
#' @example buffersteps(10, 2000, 5)
buffersteps <- function(mn=1, mx, num_step=5){
  bufs <- matrix(nrow = num_step, ncol = 1)
  v = max(1, mn)
  bufs[1] = v
  for ( i in 2:num_step) {
    v = max(v+1,(((mx-sum(bufs,na.rm=T))/v)^(1/(num_step+1 - i)) *v))
    bufs[i]=v
  }
  return (bufs)
}

#' calculates an enlarged tile bbox
#'
#' Calculates a new tile bbox that allows for buffering accross tiles and so
#' captures additional features from neighbouring tiles. Assumes earth is a
#' sphere
#'
#' @param tile.defn bbox containing columns xl, xu, yl, yu of lon/lat tile
#' @param buff buffer distance
#' @return a new tile.defn with a bbox enlarged by buffer
#'
#' @examples
#' \dontrun{
#'   xincr =0
#'   yincr =1
#'   tile.defn=data.frame(
#'   "xl"=0.5+0.041667*xincr ,
#'   "xu"=0.5+0.041667*(xincr+1),
#'   "yl"=51+0.041667*yincr ,
#'   "yu"=51+0.041667*(yincr+1)
#'   )
#'
#'   buffer_defn(tile.defn, 200)
#' }
buffer_defn <- function(tile.defn, buff=2000){
  #angle offset for buffer distance is always the same for lines of longitude
  lon_incr = buff/(2*pi*RADIUS_OF_EARTH)*360
  #limit to the coordinate space we are using - not entirely correct on a sphere but should be fine here
  yl = max(-90, tile.defn[,"yl"]-lon_incr)
  yu = min( 90, tile.defn[,"yu"]+lon_incr)

  #for latitudes we need to calculate new radii for the new longitudes
  r = radius_at_latitude(yl)
  lat_incr = buff/(2*pi*r)*360
  xl  = max(-180, tile.defn[,"xl"]-lat_incr)
  xu  = min( 180, tile.defn[,"xu"]+lat_incr)

  tile.defn["xl"] = xl
  tile.defn["xu"] = xu
  tile.defn["yl"] = yl
  tile.defn["yu"] = yu

  return(tile.defn)
}

#' find a UTM EPSG code
#'
#' Finds the EPSG code for the UTM zone that contains
#' a given lon, lat pair. Useful for converting to a
#' planar coordinate system for calculating distances etc.
#'
#' @param lon longitude coordinate in WGS84 in degrees
#' @param lat latitutde coordinate in WGS84 in degrees
#' @return EPSG code of UTM zone
lonlat2UTM = function(lon, lat) {
  utm = (floor((lon + 180) / 6) %% 60) + 1
  if(lat > 0) {
    utm + 32600
  } else{
    utm + 32700
  }
}


#' find the regular tile containing a point
#'
#' finds the bbox of the grid tile containing the given
#' lon, lat for a given regular grid division.
#'
#' @param lon longitude of contained poit
#' @param lat latitude of contained point
#' @param ref grid decomposition, on of "30 sec", "2.5 min", "15 min", "30 min", "60 min"
#' @return bbox as a vector containing xl,yl,xu, yu as elements
#'
tile_bbox <- function(lon, lat, ref="30 min") {

  part <- dplyr::case_when(
    ref == "30 sec" ~ 120,
    ref == "2.5 min" ~ 24,
    ref == "15 min" ~ 4,
    ref == "30 min" ~ 2,
    ref == "60 min" ~ 1,
    TRUE ~ 0 #trap this
  )

  xl = floor(lon*part)/part
  yl = floor(lat*part)/part
  xu = xl+1./part
  yu = yl+1./part
  return( c(xl,yl,xu,yu))
}


#' gets the centre of the given tile definition
get_center <- function(tile.defn){
  cx = tile.defn["xl"] + (tile.defn["xu"]-tile.defn["xl"])*0.5
  cy = tile.defn["yl"] + (tile.defn["yu"]-tile.defn["yl"])*0.5
  return(c(cx,cy))
}

#' converts line to well known text
#'
#' converts the data frame of coordinates representing a single line (e.g way_id)
#' into a wkt linestring
#'
#' @param sub.tbl data frame with columns lon and lat describing a single line
#' @return dataframe containing a string a linestring geometry in WKT
#'
to_linestring <- function(sub.tbl){
  x = sf::st_as_text(
    sf::st_linestring(
      sub.tbl %>% dplyr::select(lon,lat) %>% as.matrix()
    )
  )
  return (data.frame(x))
}

#' implements interface for the buffer_and_burn web service
#' parameters are encoded in json
#' @param clon (longitude of tile center)
#' @param clat (latitude of tile centre)
#' @param geom (geojson string)
#' @param bucket location of an aws bucket to store the results in
#' @return returns boolean to say whether the result was saved sucessfully
buffer_and_burn <- function(params.json){

  #----- parse the request parameters
  params.df <- jsonlite::fromJSON(params.json)
  centre_tile <- c(params.df$clon, params.df$clat)
  roads <- geojsonsf::geojson_sf(params.df$geom)
  results_bucket <- params.df$bucket #use 'bnb-output' for my aws bucket

  #----- check the tile doesn't already exist
  tbox <- geouicer::tile_bbox(centre_tile[1],centre_tile[2])
  fname <- paste(sub("\\.","_",tbox[1]), "_", sub("\\.","_",tbox[2]), ".tif", sep="")
  b<- aws.s3::get_bucket(bucket=results_bucket)
  
  if (!(result <- aws.s3::head_object(fname, b)[1])) {
  
    #calculate number of steps for incremental buffering
    steps <- geouicer::buffersteps(10,2000,5)
  
    #----- run the buffering and rasterisation
    buff_raster <- geouicer::calculate_buffer_and_burn(roads, centre_tile, steps)
  
    #----- write the output to an aws bucket

    #filename is the lower left lon lat of the tile with decimal points replaced with _ and hyphen separated
    tbox <- geouicer::tile_bbox(centre_tile[1],centre_tile[2])
    fname <- paste(sub("\\.","_",tbox[1]), "_", sub("\\.","_",tbox[2]), ".tif", sep="")
  
    result <- aws.s3::s3write_using(buff_raster, raster::writeRaster, object=fname, format="GTiff",overwrite=T, bucket=b)
  }
  
  return (result)
}


#' buffers and unions the roads and burns result to a raster
#' at the grid tile covering the lon lat contained in centre_tile.
#' @param roads collection of roads in sf format
#' @centre_tile vector containing lon, lat of the centre of the tile being processed
#' @steps vector of steps to buffer incrementally to. Can be one step.
#' @return a raster mask of buffered roads.
calculate_buffer_and_burn <- function(roads, centre_tile, buff_steps){

  #remove any lines which are invalid
  roads <- roads %>%
    as.data.frame %>%
    dplyr::filter(!is.na(sf::st_is_valid(geometry)))

  #----- buffer and rasterise

  rd_multilinestring <- sf::st_sfc(sf::st_multilinestring(roads$geometry), crs=4326)

  #convert to a multilinestring so buffers are unioned

  #Need to project to a planar coordinate system to buffer by a given distance
  #using UTM, could also have used azimuthal equidistant projection centered on centre_tile
  UTM_EPSG = geouicer::lonlat2UTM(centre_tile[1], centre_tile[2])
  rd_multilinestring = sf::st_transform(rd_multilinestring, UTM_EPSG)

  #using a chain of buffers is much faster than buffering straight to 2000m
  #where the road network is reasonably dense, need to improve strategy
  #to decide number and scaling up of strategy according to number of points and roads

  buff_mls <- sf::st_buffer(rd_multilinestring, dist=buff_steps[1])
  for ( i in 2:length(buff_steps)){
    buff_mls <- sf::st_buffer(buff_mls, dist=buff_steps[i])
  }

  #transform back to wgs84
  buff_mls <- sf::st_transform(buff_mls, 4326)

  #create a new raster to burn to, initiate with value of 1 everywhere since this will be a mask
  #takes the dimensions of the regular grid cell containing center_tile
  raster_bbox = geouicer::tile_bbox(centre_tile[1], centre_tile[2])
  buff_raster = raster::raster( nrows=60, ncols=60,
                                xmn=raster_bbox[1], ymn=raster_bbox[2],
                                xmx=raster_bbox[3], ymx=raster_bbox[4],
                                crs=4326,
                                vals = 1
  )

  #burn the buffer to a raster using mask option
  r <- raster::rasterize(as(buff_mls, "Spatial"), buff_raster, fun="max", mask=T)
  return(r)
}
