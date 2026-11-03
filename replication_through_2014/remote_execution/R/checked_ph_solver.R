# Linux-only PH adapter. Native calculations are unchanged; the accepted path
# stops before the first invalid LP step. Mac source/checkpoints remain immutable.
completion_ph_validation_policy<-function()list(
 version='checked-ph-lp-1',relative_feasibility=1e-7,bound_sum_tolerance=1e-7,
 relative_objective_gap=1e-7,minimum_slack=0,
 gap='absolute complementarity gap <= native epsilon at every accepted step',
 iteration_limit='500 accepted only when native gap criterion is satisfied',
 prefix='stop before first invalid step; no later re-entry or extrapolation')
completion_ph_native_files<-function()file.path(
 'replication_through_2014/remote_execution/validation/instrumented_ph',
 c('src/ph_diagnostic.f','src/ph_diagnostic.so','src/crqfnb_original.f',
   'src/rqfnb_original.f','build_instrumentation.py','source_hashes.json'))
completion_ph_fields<-function()c('iterations','corrector_iterations','lp_n',
 'cholesky_info','exit_phase','complementarity_gap','maxit_with_positive_gap',
 'primal_inf','primal_relative','dual_inf','dual_relative','minimum_primal_slack',
 'minimum_dual_slack','bound_sum_error','max_complementarity_product',
 'primal_objective','dual_objective','objective_gap','objective_gap_relative',
 'last_dpotrs_info','nonzero_dpotrs_count','native_epsilon','native_maxit',
 'complementarity_relative','prior_residual_nonnegative_count',
 'residual_nonnegative_count','residual_exact_zero_count','residual_near_zero_1e10_count',
 'minimum_relative_residual','residual_index_checksum','residual_absolute_sum',
 'residual_near_zero_1e8_count','nonfinite_primal_dual_state_count')
completion_ph_library<-function() {
 path<-normalizePath(completion_ph_native_files()[2L],mustWork=TRUE)
 loaded<-getLoadedDLLs()
 hits<-which(vapply(loaded,function(x)identical(normalizePath(x[['path']],mustWork=FALSE),path),logical(1)))
 if(length(hits))loaded[[hits[1L]]]else dyn.load(path)
}
completion_ph_native<-function(x,y,cen,weights,grid) {
 dll<-completion_ph_library();fields<-completion_ph_fields()
 p<-ncol(x);n<-length(y);m<-length(grid)
 x<-x*weights;y<-y*weights;cc<-as.logical(cen)
 x1<-x[cc,,drop=FALSE];y1<-y[cc];n1<-length(y1)
 z<-.Fortran('crqdiag',as.integer(n),as.integer(p),a1=as.double(t(x1)),
  c1=as.double(-y1),n1=as.integer(n1),as.double(x),as.double(y),as.double(cen),
  B=double(p*m),g=as.double(grid),m=as.integer(m),rep(1,p),double(n),
  rep(1,n),rep(1,n),wn=double(n1*9),wp=double((p+3)*p),info=integer(1),
  tr=rep(NA_real_,length(fields)*m),PACKAGE=dll[['name']],NAOK=TRUE)
 trace<-matrix(z$tr,length(fields),m,dimnames=list(fields,NULL))
 calls<-which(is.finite(trace['iterations',]))
 J<-max(0L,min(m,z$m-1L));B<--matrix(z$B,p,m)[,seq_len(J),drop=FALSE]
 rownames(B)<-colnames(x)
 list(beta=B,tau=grid[seq_len(J)],native_info=z$info,native_m=z$m,
   trace=trace[,calls,drop=FALSE],trace_grid_index=calls,trace_tau=grid[calls])
}
completion_ph_validate_native<-function(z,policy=completion_ph_validation_policy()) {
 stopifnot(identical(policy,completion_ph_validation_policy()),is.matrix(z$trace),
   identical(rownames(z$trace),completion_ph_fields()),
   ncol(z$trace)==length(z$trace_grid_index),length(z$trace_tau)==ncol(z$trace),
   identical(z$trace_grid_index,seq_len(ncol(z$trace))),ncol(z$beta)==length(z$tau),
   ncol(z$trace)>=ncol(z$beta))
 t<-z$trace;reason<-rep('valid',ncol(t))
 required<-c('iterations','corrector_iterations','lp_n','cholesky_info','exit_phase',
   'complementarity_gap','maxit_with_positive_gap','primal_inf','primal_relative',
   'dual_inf','dual_relative','minimum_primal_slack','minimum_dual_slack',
   'bound_sum_error','max_complementarity_product','primal_objective','dual_objective',
   'objective_gap','objective_gap_relative','last_dpotrs_info','nonzero_dpotrs_count',
   'native_epsilon','native_maxit','complementarity_relative','nonfinite_primal_dual_state_count')
 for(i in seq_len(ncol(t))) {
   r<-t[,i]
   reason[i]<-if(!is.finite(r['cholesky_info'])||!is.finite(r['exit_phase']))'invalid_native_status'else
     if(r['cholesky_info']!=0||r['exit_phase']%in%c(0,1))'native_cholesky_failure'else
     if(r['exit_phase']!=2)'invalid_exit_phase'else
     if(any(!is.finite(r[required]))||r['nonfinite_primal_dual_state_count']!=0)'nonfinite_lp_state'else
     if(r['last_dpotrs_info']!=0||r['nonzero_dpotrs_count']!=0)'native_dpotrs_failure'else
     if(r['native_epsilon']!=1e-8||r['native_maxit']!=500||r['iterations']<0||
        r['iterations']>r['native_maxit']||r['iterations']!=floor(r['iterations']))'invalid_native_controls'else
     if(r['maxit_with_positive_gap']!=0||
        (r['iterations']>=r['native_maxit']&&abs(r['complementarity_gap'])>r['native_epsilon']))'iteration_limit_unconverged'else
     if(abs(r['complementarity_gap'])>r['native_epsilon'])'gap_not_converged'else
     if(r['minimum_primal_slack']<policy$minimum_slack||r['minimum_dual_slack']<policy$minimum_slack)'negative_lp_slack'else
     if(r['bound_sum_error']<0||r['bound_sum_error']>policy$bound_sum_tolerance)'invalid_bound_sum'else
     if(r['primal_relative']<0||r['dual_relative']<0||
        r['primal_relative']>policy$relative_feasibility||r['dual_relative']>policy$relative_feasibility)'lp_feasibility_failure'else
     if(r['objective_gap_relative']<0||r['objective_gap_relative']>policy$relative_objective_gap)'lp_objective_gap_failure'else'valid'
 }
 # Trace records are necessary even when the returned coefficient is finite.
 bad<-which(reason!='valid')[1L]
 accepted<-min(ncol(z$beta),if(is.na(bad))ncol(z$beta)else bad-1L)
 list(reason=reason,first_invalid_index=if(is.na(bad))NA_integer_ else as.integer(bad),
   first_invalid_tau=if(is.na(bad))NA_real_ else z$trace_tau[bad],
   first_invalid_reason=if(is.na(bad))'none'else reason[bad],
   accepted_steps=accepted,returned_steps=ncol(z$beta),truncated=accepted<ncol(z$beta),policy=policy)
}
completion_checked_ph_fit<-function(prepared,weights,coefficient_names=NULL,
                                    weight_constant=NULL,keep_path=FALSE) {
 state<-new.env(parent=emptyenv())
 checked_native<-function(x,y,cen,weights,grid) {
   z<-completion_ph_native(x,y,cen,weights,grid)
   checked<-completion_ph_validate_native(z)
   state$raw<-z;state$check<-checked
   take<-seq_len(checked$accepted_steps)
   list(beta=z$beta[,take,drop=FALSE],tau=z$tau[take],
     native_info=z$native_info,native_m=z$native_m)
 }
 # Copy only the function's environment, preserving its exact calculations and
 # original input/rank checks. The globally frozen native function is untouched.
 fn<-revised_cqr_fit;private<-new.env(parent=environment(fn))
 private$revised_pen_native<-checked_native;environment(fn)<-private
 ans<-fn(prepared,weights,coefficient_names,weight_constant,keep_path)
 ans$diagnostics$solver_version<-'checked-pen-lp-1.0'
 ans$diagnostics$lp_validation_policy<-completion_ph_validation_policy()
 if(exists('check',envir=state,inherits=FALSE)) {
   ck<-state$check;tr<-state$raw$trace
   ans$diagnostics$lp_steps_checked<-ncol(tr)
   ans$diagnostics$lp_steps_accepted<-ck$accepted_steps
   ans$diagnostics$lp_returned_steps<-ck$returned_steps
   ans$diagnostics$lp_first_invalid_index<-ck$first_invalid_index
   ans$diagnostics$lp_first_invalid_tau<-ck$first_invalid_tau
   ans$diagnostics$lp_first_invalid_reason<-ck$first_invalid_reason
   ans$diagnostics$lp_path_truncated<-ck$truncated
   ans$diagnostics$lp_invalid_steps<-sum(ck$reason!='valid')
   ans$diagnostics$lp_first_invalid_step<-if(is.na(ck$first_invalid_index))NULL else
     list(grid_index=ck$first_invalid_index,tau=ck$first_invalid_tau,
       reason=ck$first_invalid_reason,metrics=tr[,ck$first_invalid_index])
   ans$diagnostics$lp_accepted_tau_max<-if(ck$accepted_steps)state$raw$tau[ck$accepted_steps]else NA_real_
   accepted<-seq_len(ck$accepted_steps)
   metric<-function(field)if(length(accepted))max(abs(tr[field,accepted]))else NA_real_
   ans$diagnostics$lp_max_iterations<-metric('iterations')
   ans$diagnostics$lp_max_absolute_gap<-metric('complementarity_gap')
   ans$diagnostics$lp_max_primal_relative<-metric('primal_relative')
   ans$diagnostics$lp_max_dual_relative<-metric('dual_relative')
   ans$diagnostics$lp_max_objective_gap_relative<-metric('objective_gap_relative')
   if(ck$truncated)ans$diagnostics$failure_class<-paste0('checked_',ck$first_invalid_reason)
   if(keep_path) {
     ans$validation_trace<-data.frame(grid_index=state$raw$trace_grid_index,
       tau=state$raw$trace_tau,t(tr),validation_reason=ck$reason,check.names=FALSE)
     if(ck$truncated)ans$unchecked_native_path<-list(tau=state$raw$tau,
       coefficients=revised_cqr_backtransform(prepared,state$raw$beta),
       native_info=state$raw$native_info)
   }
 }else {
   ans$diagnostics$lp_steps_checked<-0L;ans$diagnostics$lp_steps_accepted<-0L
   ans$diagnostics$lp_first_invalid_reason<-'not_called_input_or_native_error'
   ans$diagnostics$lp_path_truncated<-FALSE
 }
 ans
}
