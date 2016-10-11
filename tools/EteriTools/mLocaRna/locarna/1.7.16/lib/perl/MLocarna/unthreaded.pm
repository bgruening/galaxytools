## just do the same as foreach_par but sequentially
##
sub foreach_par {
    my ($sub_ref,$argument_lists_ref,$thread_num) = @_;
    
    my @argument_lists = @{ $argument_lists_ref };
    
    for (my $my_job=0; $my_job <= $#argument_lists ; $my_job++) { 
	## start job
	$sub_ref->(@{ $argument_lists[$my_job] });
    }
}


## generate a thread-unique filename
sub threadsafe_name($) {
    my ($name) = @_;
    return $name;
}


#do nothing
sub share(\[$@%]) {
    return;
}

return 1;
