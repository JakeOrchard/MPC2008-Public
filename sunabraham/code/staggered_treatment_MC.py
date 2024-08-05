import pandas as pd
import numpy as np

def staggered_treatment_MC(T=2, TC=2, Nset=[10,10], treatvals=[1,1], treattime=[0,1], 
                           sigmaT = 0.01, sigmai = 0.01, sigmaiT = 0.01, sigmaC = 0.01):


    # convert to arrays
    if len(Nset)<len(treattime):
        Nset = Nset * np.ones(len(treattime))
    else:
        Nset = np.array(Nset)

    treatvals = np.array(treatvals)
    treattime = np.array(treattime)   

    data = np.array([], dtype=int)

    # Date and TFE are constant across cohorts
    Date = np.array(range(T)).reshape([T,1])     
    TFE = np.random.normal(size=[T,1]) * sigmaT

    # constructing cohorts here
    for tstart in range(1+ T-TC):

        CohortDate = Date[tstart:tstart+TC]

        CohortTFE = TFE[tstart:tstart+TC]

        for i, N in enumerate(Nset):

            CFE = np.random.normal(size=[TC,1]) * sigmaC

            D = np.zeros([TC, 1])

            if i<=len(treattime)-1:

                D[treattime[i],] = 1

                treatval = treatvals[i]
                                 
            else:
                treatval = 0


            for n in range(N):

                HHind = np.ones([TC,1]) * (n + Nset[0:i].sum() + tstart*Nset.sum())

                c = np.random.normal() * sigmai

                Y = c + CFE + CohortTFE + treatval * D + np.random.normal(size=[TC,1]) * sigmaiT

                datatoadd = np.concatenate((CohortDate, HHind, Y, D), axis=1)

                data = np.concatenate((data, datatoadd)) if data.size else datatoadd

    df = pd.DataFrame(data, 
                      columns=['Date', 'HHind', 'Y', 'D'])

    df.Date = df.Date.astype(np.int64)
    df.HHind = df.HHind.astype(np.int64)

    if df.index.is_unique==False:
        print('Warning: Index not unique.')

    return df