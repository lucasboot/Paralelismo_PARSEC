double pgain ( long x, Points *points, double z, long int *numcenters )
{
    int i;
    int number_of_centers_to_close = 0;

    static double *work_mem;
    static double gl_cost_of_opening_x;
    static int gl_number_of_centers_to_close;

    int stride = *numcenters + 2;
    //make stride a multiple of CACHE_LINE
    int cl = CACHE_LINE/sizeof ( double );
    if ( stride % cl != 0 ) {
        stride = cl * ( stride / cl + 1 );
    }
    int K = stride - 2 ; // K==*numcenters

    //my own cost of opening x
    double cost_of_opening_x = 0;
    work_mem = ( double* ) malloc ( 2 * stride * sizeof ( double ) );
    gl_cost_of_opening_x = 0;
    gl_number_of_centers_to_close = 0;

    int count = 0;
    //my *lower* fields
    double* lower;
    //global *lower* fields
    double* gl_lower;

    #pragma omp parallel
    {
        /*
         * For each center, we have a *lower* field that indicates
         * how much we will save by closing the center.
         */

        int i;
        #pragma omp for private(i)
        for ( i = 0; i < points->num; i++ ) {
            if ( is_center[i] ) {
                #pragma omp critical
                center_table[i] = count++;
            }
        }

        #pragma omp single
        work_mem[0] = 0;

        #pragma omp sections
        {
            //now we finish building the table. clear the working memory.
            #pragma omp section
            memset ( switch_membership, 0, points->num * sizeof ( bool ) );

            #pragma omp section
            memset ( work_mem, 0, stride*sizeof ( double ) );

            #pragma omp section
            memset ( work_mem+stride,0,stride*sizeof ( double ) );
        }

        #pragma omp single
        {
            lower = &work_mem[0];
            gl_lower = &work_mem[stride];
        }

        float x_cost, current_cost;

        #pragma omp for private(i, x_cost, current_cost)
        for ( i = 0; i < points->num; i++ ) {

            x_cost = dist ( points->p[i], points->p[x], points->dim ) * points->p[i].weight;
            current_cost = points->p[i].cost;

            if ( x_cost < current_cost ) {

                // point i would save cost just by switching to x
                // (note that i cannot be a median,
                // or else dist(p[i], p[x]) would be 0)

                switch_membership[i] = 1;

                #pragma omp critical
                cost_of_opening_x += x_cost - current_cost;

            } else {

                // cost of assigning i to x is at least current assignment cost of i

                // consider the savings that i's **current** median would realize
                // if we reassigned that median and all its members to x;
                // note we've already accounted for the fact that the median
                // would save z by closing; now we have to subtract from the savings
                // the extra cost of reassigning that median and its members
                int assign = points->p[i].assign;

                #pragma omp critical
                lower[center_table[assign]] += current_cost - x_cost;
            }
        }

        // at this time, we can calculate the cost of opening a center
        // at x; if it is negative, we'll go through with opening it

        double low;

        #pragma omp for private(i, low)
        for ( int i = 0; i < points->num; i++ ) {
            if ( is_center[i] ) {
                low = z + work_mem[center_table[i]];
                #pragma omp critical
                gl_lower[center_table[i]] = low;
                if ( low > 0 ) {
                    // i is a median, and
                    // if we were to open x (which we still may not) we'd close i

                    // note, we'll ignore the following quantity unless we do open x
                    #pragma omp atomic
                    ++number_of_centers_to_close;
                    #pragma omp critical
                    cost_of_opening_x -= low;
                }
            }
        }

        #pragma omp sections
        {
            //use the rest of working memory to store the following
            #pragma omp section
            work_mem[K] = number_of_centers_to_close;

            #pragma omp section
            work_mem[K+1] = cost_of_opening_x;

            #pragma omp section
            gl_number_of_centers_to_close = ( int ) work_mem[K];

            #pragma omp section
            gl_cost_of_opening_x = z + work_mem[K+1];
        }

        // Now, check whether opening x would save cost; if so, do it, and
        // otherwise do nothing
        bool close_center;

        if ( gl_cost_of_opening_x < 0 ) {
            //  we'd save money by opening x; we'll do it
            #pragma omp for private(i)
            for ( i = 0; i < points->num; i++ ) {
                close_center = gl_lower[center_table[points->p[i].assign]] > 0 ;
                if ( switch_membership[i] || close_center ) {
                    // Either i's median (which may be i itself) is closing,
                    // or i is closer to x than to its current median
                    points->p[i].cost = points->p[i].weight * dist ( points->p[i], points->p[x], points->dim );
                    points->p[i].assign = x;
                }
            }
            #pragma omp for private(i)
            for ( i = 0; i < points->num; i++ ) {
                if ( is_center[i] && gl_lower[center_table[i]] > 0 ) {
                    is_center[i] = false;
                }
            }
            if ( x >= 0 && x < points->num ) {
                is_center[x] = true;
            }

            #pragma omp single
            *numcenters = *numcenters + 1 - gl_number_of_centers_to_close;
        } else {
            #pragma omp single
            gl_cost_of_opening_x = 0;  // the value we'll return
        }

        #pragma omp single
        free ( work_mem );
    }

    return -gl_cost_of_opening_x;
}