library(geouicer)

testfile = "~/Documents/GitHub/athena-buildings/output/test.geojson"
gjson = file(testfile)
roads <- geojsonsf::geojson_sf(gjson)
centre_tile = c( 1.16375, 51.1145)
bsteps = geouicer::buffersteps(mn=15,mx=2000,num_step=5)
r = geouicer::calculate_buffer_and_burn(roads, centre_tile, bsteps)

#check the bbox using xmin etc.
plot(r)




params.df <- data.frame(clon = centre_tile[1], 
                        clat = centre_tile[2], 
                        geom=testfile,
                        bucket='bnb-output')

params.json <- jsonlite::toJSON(params.df)

buffer_and_burn(params.json)
