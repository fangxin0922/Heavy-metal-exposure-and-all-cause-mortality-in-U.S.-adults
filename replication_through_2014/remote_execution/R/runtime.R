# Explicit runtime provenance for the isolated Linux execution. Scientific
# prepared artifacts and the original macOS run remain immutable evidence.
remote_library_path<-Sys.getenv('HM_REMOTE_R_LIBRARY',file.path(normalizePath('.'),'replication_10_24/library'))
optional_library_path<-Sys.getenv('HM_R_LIBRARY',unset='')
.libPaths(c(remote_library_path,optional_library_path[nzchar(optional_library_path)],.libPaths()))
remote_environment_ready<-function() {
 stopifnot(file.exists('replication_through_2014/remote_execution/validation/environment_ready.txt'))
 invisible(TRUE)
}
remote_runtime_record<-function(require_linux=FALSE) {
 info<-sessionInfo()
 blas<-normalizePath(info$BLAS,mustWork=TRUE)
 lapack<-normalizePath(La_library(),mustWork=TRUE)
 if(require_linux)stopifnot(identical(Sys.info()[['sysname']],'Linux'),
   grepl('openblas',blas,ignore.case=TRUE),
   !nzchar(Sys.getenv('LD_PRELOAD')),
   Sys.getenv('OPENBLAS_NUM_THREADS')=='1',Sys.getenv('OMP_NUM_THREADS')=='1')
 packages<-c('quantreg','survival','SparseM','MatrixModels','Matrix','stats','methods')
 pkg_files<-unlist(lapply(packages,function(pkg) {
   p<-find.package(pkg)
   c(file.path(p,'DESCRIPTION'),list.files(file.path(p,'R'),pattern='[.](rdb|rdx)$',full.names=TRUE),
     list.files(file.path(p,'libs'),recursive=TRUE,full.names=TRUE))
 }))
 runtime_files<-c(blas,lapack,file.path(R.home('bin'),'exec','R'),
   list.files(R.home('lib'),pattern='^libR([.]|blas|lapack)',full.names=TRUE),pkg_files)
 checked_ph_files<-file.path('replication_through_2014/remote_execution/validation/instrumented_ph/src',
   c('ph_diagnostic.f','ph_diagnostic.so'))
 if(require_linux)stopifnot(all(file.exists(checked_ph_files)))
 checked_ph_files<-checked_ph_files[file.exists(checked_ph_files)]
 runtime_files<-c(runtime_files,checked_ph_files)
 runtime_files<-runtime_files[!dir.exists(runtime_files)]
 runtime_files<-unique(normalizePath(runtime_files,mustWork=TRUE))
 stopifnot(length(runtime_files)>0L,all(file.exists(runtime_files)))
 list(R_version=R.version.string,R_platform=R.version$platform,
   os=unname(Sys.info()[c('sysname','release','machine')]),BLAS=blas,LAPACK=lapack,
   package_versions=setNames(vapply(packages,function(p)as.character(packageVersion(p)),character(1)),packages),
   checked_ph_native_hashes=tools::md5sum(normalizePath(checked_ph_files,mustWork=TRUE)),
   thread_environment=Sys.getenv(c('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','OPENBLAS_CORETYPE','LD_PRELOAD')),
   hashes=tools::md5sum(runtime_files))
}
