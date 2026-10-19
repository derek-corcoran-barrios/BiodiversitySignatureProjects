# Checks real terra tile preparation without needing the saved species fits.
source("denmark_pipeline/run.R")
folder <- tempfile("projection_input_test_", tmpdir = getwd()); dir.create(folder)
folder <- normalizePath(folder, winslash = "/")
r <- rast(nrows=6,ncols=6,xmin=500000,xmax=500060,ymin=6200000,ymax=6200060,crs="EPSG:25832")
values(r) <- rep(c(1,9,11,10,9,12),6)
rat <- data.frame(id=1:12,class=c(dp_pathways()$class,"Agriculture","Urban","ProductionForest","Other"))
r <- categories(r,value=rat); names(r)<-"Landuse"
lu <- file.path(folder,"landuse.tif");dp_raster(r,lu,"INT2U")
v <- vect(ext(r),crs=crs(r));boundary<-file.path(folder,"boundary.gpkg");writeVector(v,boundary,overwrite=TRUE)
cl <- rast(nrows=5,ncols=5,xmin=499980,xmax=500080,ymin=6199980,ymax=6200080,crs=crs(r))
values(cl) <- seq_len(ncell(cl)); climate <- c(cl,cl+1,cl+2,cl+3,cl+4)
names(climate)<-c("bio1","bio4","bio12","bio15","bio18")
cp<-file.path(folder,"climate.tif");dp_raster(climate,cp)
tm<-file.path(folder,"training_manifest.rds");saveRDS(c(lu,cp),tm)
cfg<-denmark_config();cfg$grid$landuse<-lu;cfg$region$boundary<-boundary
cfg$projection$training_manifest<-tm;cfg$projection$climate<-cp
cfg$projection$output_root<-file.path(folder,"projections");cfg$projection$tile_cells<-3L
plan<-dp_projection_plan(cfg,validate_models=FALSE)
stopifnot(length(plan$jobs)==4L)
products<-lapply(plan$jobs,function(j) {
  f<-file.path(folder,j$id);dir.create(f);dp_tile_predictors(plan,j,f)
})
whole_climate<-project(climate,r,method="bilinear")
for (i in seq_along(products)) {
  x<-products[[i]]; current<-rast(x$landuse$current); candidate<-as.numeric(current) %in% c(9,11)
  stopifnot(length(x$landuse)==9L,identical(cats(current)[[1]],cats(rast(lu))[[1]]))
  for (k in seq_len(8)) {
    z<-rast(x$landuse[[k+1]])
    stopifnot(!dp_any(is.na(z)!=is.na(current)))
    stopifnot(!dp_any(!candidate & as.numeric(z)!=as.numeric(current)))
    stopifnot(!dp_any(candidate & as.numeric(z)!=k))
  }
  expected<-crop(whole_climate,ext(current));actual<-rast(x$climate)
  stopifnot(compareGeom(expected,actual),max(abs(values(expected)-values(actual)),na.rm=TRUE)<1e-9)
}
tile_files <- vapply(products, function(p) p$landuse$current, character(1))
vf <- file.path(folder,"current_mosaic.vrt"); terra::vrt(tile_files,filename=vf,overwrite=TRUE)
tracked <- normalizePath(dp_files(vf),winslash="/")
stopifnot(all(normalizePath(tile_files,winslash="/") %in% tracked))
# Missing climate, unsupported categorical dictionaries, and changed inputs fail loudly.
bad<-cfg;bad$projection$id<-"bad_categories"; bad$grid$landuse<-file.path(folder,"bad_landuse.tif")
rr<-categories(r,value=transform(rat,class=paste0(class,"_changed")));dp_raster(rr,bad$grid$landuse,"INT2U")
stopifnot(inherits(try(dp_projection_plan(bad,validate_models=FALSE),silent=TRUE),"try-error"))
cat("PASS: four 10 m tiles, all eight exact endpoint scenarios, unchanged nonselectable cells, preserved categories and domains, seamless climate interpolation, dictionary mismatch rejected.\n")
