# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 14:13:02 2020

@author: eliza

Master document for the ODEs of a co-culture system where one cell
expressed heterologous protein ea, the other expresses heterologous protein eb,
and both are required for starch degradation in the extracellular environment
"""


def dydt_duo(t,y,p, *args, **kwds):
    cell_a_size= 21
    cell_b_size= 21
    y_a= y[:cell_a_size]
    y_b= y[cell_a_size:(cell_a_size+cell_b_size)]
    y_ext= y[(cell_a_size+cell_b_size):]
    dydt_a = dydt_duo_ea(t,y,p,y_a,y_ext) 
    dydt_b = dydt_duo_eb(t,y,p,y_b,y_ext)
    dydt_ext = dydt_duo_ext(t,y,p,y_a,y_b,y_ext, *args, **kwds)
    return(dydt_a + dydt_b + dydt_ext)


#CELL A

def dydt_duo_ea (t, y, p, y_a, y_ext):
    
    ma = y_a[0]
    rma = y_a[1]
    ea = y_a[2] 
    
    mq= y_a[3]
    rmq= y_a[4]
    q= y_a[5]
    
    mm = y_a[6]
    rmm = y_a[7]
    em = y_a[8]
    
    mt = y_a[9]
    rmt= y_a[10]
    et= y_a[11]
    
    mr = y_a[12]
    rmr = y_a[13]
    r = y_a[14]
    
    si= y_a[15]
    a= y_a[16]
    
    mc = y_a[17]
    rmc = y_a[18]
    ec = y_a[19] 
    
    z = y_a[20]
    
    #from extracellular ODEs
    s= y_ext[2]
    
    #variable constants
    gamma= (p["gmax"]*a/(p["Kgamma"] + a))
    ttrate = ((rmq + rmr + rmt + rmm + rma + rmc)*gamma)
    lam = ttrate/p["M"]
    rate_import = p["v_et"]*et*s/(p["Km_et"] + s)
    rate_metabolism = p["v_em"]*em*si/(p["Km_em"] + si)
    
    #define ODEs
    
    #enzyme a
    
    #ma enzyme a mRNA
    dmadt= (p["w_ea"]*a/(p["theta_x"] + a)) + p["ku"]*rma + gamma/p["n_a"]*rma - p["kb"]*r*ma - p["dm"]*ma - lam*ma
    #rma ribo:mRNA complex 
    drmadt= p["kb"]*r*ma - p["ku"]*rma - gamma/p["n_a"]*rma - lam*rma
    #enzyme alpha (a-amylase)
    deadt= gamma/p["n_a"]*rma - lam*ea - p["d_ea"]*ea
    
    #housekeeping
    
    #mq housekeeping protein mRNA
    dmqdt= (p["w_q"]*a/(p["theta_x"] + a)) * (1/(1 + (q/p["Kq"])**p["nq"])) + (p["ku"]*rmq) + (rmq*(gamma/p["n_x"])) - (p["kb"]*r*mq) - (p["dm"]*mq) - (lam*mq)
    #rmq ribo:housekeeping mRNA
    drmqdt= (p["kb"]*r*mq) - (p["ku"]*rmq) - (rmq*(gamma/p["n_x"])) - (lam*rmq)
    #q housekeeping protein
    dqdt= (rmq*(gamma/p["n_x"])) - (lam*q)
    
    #metabolic
    
    #mm metabolic enzyme mRNA
    dmmdt= (p["w_m"]*a/(p["theta_x"] + a)) + p["ku"]*rmm + gamma/p["n_x"]*rmm - p["kb"]*r*mm - p["dm"]*mm - lam*mm
    #rmm Ribo:metabolic enzyme mRNA
    drmmdt= p["kb"]*r*mm - p["ku"]*rmm - gamma/p["n_x"]*rmm - lam*rmm
    #em enzyme metabolic
    demdt= gamma/p["n_x"]*rmm - lam*em
    
    #transport
    
    #mt transport enzyme mRNA
    dmtdt= (p["w_t"]*a/(p["theta_x"] + a)) + p["ku"]*rmt + gamma/p["n_x"]*rmt - p["kb"]*r*mt - p["dm"]*mt - lam*mt
    #rmt Ribo:transport enzyme mRNA
    drmtdt= p["kb"]*r*mt - p["ku"]*rmt - gamma/p["n_x"]*rmt - lam*rmt    
    #et enzyme transport
    detdt= gamma/p["n_x"]*rmt - lam*et

    #ribosomes

    #mr ribo mRNA
    dmrdt= (p["w_r"]*a/(p["theta_r"] + a)) - (p["kb"]*r*mr) + (p["ku"]*rmr) + rmr*(gamma/p["n_r"]) - (p["dm"]*mr) - (lam*mr)
    #rmr ribo:ribo mRNA
    drmrdt= (p["kb"]*r*mr) - (p["ku"]*rmr) - rmr*(gamma/p["n_r"]) - (lam*rmr)
    #r ribosomes
             #ribsome production     #ribosome unbinding of all mrna complexes                             #completion of translation for all proteins                                                                         #binding of free ribosome to mrna                                        #dilution
    drdt=  (rmr*(gamma/p["n_r"])) + p["ku"]*rma + p["ku"]*rmc + p["ku"]*rmq + p["ku"]*rmm + p["ku"]*rmt + p["ku"]*rmr + rma*(gamma/p["n_a"]) + rmc*(gamma/p["n_c"]) + rmq*(gamma/p["n_x"]) + rmm*(gamma/p["n_x"]) + rmt*(gamma/p["n_x"]) + rmr*(gamma/p["n_r"]) - p["kb"]*r*ma - p["kb"]*r*mc - p["kb"]*r*mq - p["kb"]*r*mm - p["kb"]*r*mt - p["kb"]*r*mr - lam*r
    
    #intracellular energy
    
    #Imported glucose 
    dsidt= rate_import - rate_metabolism - lam*si 
    #a available energy
    dadt= p["ns"]*rate_metabolism - ttrate - lam*a
    
    #mx enzyme c mRNA
    dmcdt= (p["w_ec"]*a/(p["theta_x"] + a)) + p["ku"]*rmc + gamma/p["n_x"]*rmc - p["kb"]*r*mc - p["dm"]*mc - lam*mc
    #rmc ribo:mRNA complex 
    drmcdt= p["kb"]*r*mc - p["ku"]*rmc - (gamma/p["n_x"])*rmc - lam*rmc
    #enzyme c (arbitrary enzyme for desirable metabolic pathway)
    decdt= (gamma/p["n_c"])*rmc - lam*ec - p["dg_ec"]*ec    

    #reaction converting "energy" (ATP) into a product
    
    dzdt= p["v_ec"] * ec * a / (p["Km_ec"] + a) - lam*z

    return (dmadt, drmadt, deadt, dmqdt, drmqdt, dqdt, dmmdt, drmmdt, demdt, dmtdt, drmtdt, detdt, dmrdt, drmrdt, drdt, dsidt, dadt, dmcdt, drmcdt, decdt, dzdt)


#CELL B

def dydt_duo_eb (t, y, p, y_b, y_ext):
    
    mb = y_b[0] 
    rmb= y_b[1]
    eb= y_b[2]
    
    mq= y_b[3]
    rmq= y_b[4]
    q= y_b[5]
    
    mm = y_b[6]
    rmm = y_b[7]
    em = y_b[8]
    
    mt = y_b[9]
    rmt= y_b[10]
    et= y_b[11]
    
    mr = y_b[12]
    rmr = y_b[13]
    r = y_b[14]
    
    si= y_b[15]
    a= y_b[16]
    
    mc = y_b[17] 
    rmc= y_b[18]
    ec = y_b[19]
    
    z = y_b[20]

    #from extracellular ODEs
    s= y_ext[2]
    
    #variable constants
    gamma= (p["gmax"]*a/(p["Kgamma"] + a))
    ttrate = ((rmq + rmr + rmt + rmm + rmb)*gamma)
    lam = ttrate/p["M"]
    rate_import = p["v_et"]*et*s/(p["Km_et"] + s)
    rate_metabolism = p["v_em"]*em*si/(p["Km_em"] + si)
    
    #define ODEs
    
    #enzyme b
    
    #mb enzyme b mRNA
    dmbdt= (p["w_eb"]*a/(p["theta_x"] + a)) + p["ku"]*rmb + gamma/p["n_b"] * rmb - p["kb"]*r*mb - p["dm"]*mb - lam*mb
    #rmb ribo:mRNA complex
    drmbdt= p["kb"]*r*mb - p["ku"]*rmb - gamma/p["n_b"]*rmb - lam*rmb
    #enzyme b (glucoamylase)
    debdt= gamma/p["n_b"]*rmb - lam*eb - p["d_eb"]*eb
    
    #housekeeping
    
    #mq housekeeping protein mRNA
    dmqdt= (p["w_q"]*a/(p["theta_x"] + a)) * (1/(1 + (q/p["Kq"])**p["nq"])) + (p["ku"]*rmq) + (rmq*(gamma/p["n_x"])) - (p["kb"]*r*mq) - (p["dm"]*mq) - (lam*mq)
    #rmq ribo:housekeeping mRNA
    drmqdt= (p["kb"]*r*mq) - (p["ku"]*rmq) - (rmq*(gamma/p["n_x"])) - (lam*rmq)
    #q housekeeping protein
    dqdt= (rmq*(gamma/p["n_x"])) - (lam*q)
    
    #metabolic
    
    #mm metabolic enzyme mRNA
    dmmdt= (p["w_m"]*a/(p["theta_x"] + a)) + p["ku"]*rmm + gamma/p["n_x"]*rmm - p["kb"]*r*mm - p["dm"]*mm - lam*mm
    #rmm Ribo:metabolic enzyme mRNA
    drmmdt= p["kb"]*r*mm - p["ku"]*rmm - gamma/p["n_x"]*rmm - lam*rmm
    #em enzyme metabolic
    demdt= gamma/p["n_x"]*rmm - lam*em
    
    #transport
    
    #mt transport enzyme mRNA
    dmtdt= (p["w_t"]*a/(p["theta_x"] + a)) + p["ku"]*rmt + gamma/p["n_x"]*rmt - p["kb"]*r*mt - p["dm"]*mt - lam*mt
    #rmt Ribo:transport enzyme mRNA
    drmtdt= p["kb"]*r*mt - p["ku"]*rmt - gamma/p["n_x"]*rmt - lam*rmt    
    #et enzyme transport
    detdt= gamma/p["n_x"]*rmt - lam*et

    #ribosomes

    #mr ribo mRNA
    dmrdt= (p["w_r"]*a/(p["theta_r"] + a)) - (p["kb"]*r*mr) + (p["ku"]*rmr) + rmr*(gamma/p["n_r"]) - (p["dm"]*mr) - (lam*mr)
    #rmr ribo:ribo mRNA
    drmrdt= (p["kb"]*r*mr) - (p["ku"]*rmr) - rmr*(gamma/p["n_r"]) - (lam*rmr)
    #r ribosomes
                 #ribsome production     #ribosome unbinding of all mrna complexes                             #completion of translation for all proteins                                                                         #binding of free ribosome to mrna                                        #dilution
    drdt=  (rmr*(gamma/p["n_r"])) + p["ku"]*rmb + p["ku"]*rmc + p["ku"]*rmq + p["ku"]*rmm + p["ku"]*rmt + p["ku"]*rmr + rmb*(gamma/p["n_b"]) + rmc*(gamma/p["n_c"]) + rmq*(gamma/p["n_x"]) + rmm*(gamma/p["n_x"]) + rmt*(gamma/p["n_x"]) + rmr*(gamma/p["n_r"]) - p["kb"]*r*mb - p["kb"]*r*mc - p["kb"]*r*mq - p["kb"]*r*mm - p["kb"]*r*mt - p["kb"]*r*mr - lam*r
    
    #intracellular energy
    
    #Imported glucose 
    dsidt= rate_import - rate_metabolism - lam*si 
    #a available energy
    dadt= p["ns"]*rate_metabolism - ttrate - lam*a
    
    #mx enzyme c mRNA
    dmcdt= (p["w_ec"]*a/(p["theta_x"] + a)) + p["ku"]*rmc + gamma/p["n_x"]*rmc - p["kb"]*r*mc - p["dm"]*mc - lam*mc
    #rmc ribo:mRNA complex 
    drmcdt= p["kb"]*r*mc - p["ku"]*rmc - (gamma/p["n_x"])*rmc - lam*rmc
    #enzyme c (arbitrary enzyme for desirable metabolic pathway)
    decdt= (gamma/p["n_c"])*rmc - lam*ec - p["dg_ec"]*ec
    
    #reaction converting "energy" (ATP) into a product
    
    dzdt= p["v_ec"] * ec * a / (p["Km_ec"] + a) - lam*z


    return (dmbdt, drmbdt, debdt, dmqdt, drmqdt, dqdt, dmmdt, drmmdt, demdt, dmtdt, drmtdt, detdt, dmrdt, drmrdt, drdt, dsidt, dadt, dmcdt, drmcdt, decdt, dzdt)


#EXTRACELLULAR

def dydt_duo_ext (t, y, p, y_a, y_b, y_ext, substrate):
    
    s0 = y_ext[0]
    s1 = y_ext[1]
    s = y_ext[2]
    
    #from intracellular ODEs
    ea= y_a[2]
    eb= y_b[2]
    et_a = y_a[11] 
    et_b = y_b[11]
    
    rate_import_a = p["v_et"]*et_a*s/(p["Km_et"] + s)
    rate_import_b = p["v_et"]*et_b*s/(p["Km_et"] + s)
    
    #if glucose is the primary substrate: 
    if substrate == "glucose": 
        ds0dt = 0
        ds1dt = 0      
        dsdt = p["s_in"] - (rate_import_a * p["Na"]) - (rate_import_b * p["Nb"]) - p["d_s"]*s
        
    #if glucose is not the primary substrate then the primary substrate is
    #starch and the rates will change
    elif substrate == "starch":
        
        ###determine the extent of degradation so far - using Fujii et al equations
        eta = ( s / (p["xi"] * s0) ) * ( p["mw_glc"] / p["mw_s0"] )
       
        
        #depending on extent of degradation the rates of glucoamylase will change
        if eta < 0.4:
            #1st degradation - action of a-amylase
            rate_ea = p["v_ea"] * ea * s0 / (p["Km_ea"] + s0) 
            #2nd degradation - action of glucoamylase
            rate_eb = p["v_eb"]  * eb * s1 / (p["Km_eb"] * (1 + s / p["K_i"]) + s1)
    
        else:
            print("eta more than 0.4")
            #rate/Km of glucoamylase changes with size of molecules
            v_m2 = p["v_eb"] * eb * (1 - 0.28*(eta - 0.4))
            K_m2 = p["Km_eb"] * (1 + 0.69*(eta - 0.4))
            
            rate_ea = 0 
            rate_eb = v_m2 * s1 / (K_m2 * (1 + s / p["K_i"]) + s1)
    

        ds0dt = p["s0_in"] - rate_ea * p["Na"] - p["d_s0"]*s0
        ds1dt= rate_ea * p["Na"] - rate_eb * p["Nb"] - p["d_s1"]*s1
        dsdt = rate_eb * p["Nb"] - rate_import_a * p["Na"] - rate_import_b * p["Nb"] - p["d_s"]*s
        
    return (ds0dt, ds1dt, dsdt)
