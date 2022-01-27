# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 16:59:26 2020

@author: Eliza

ODEs for a monoculture system expressing 2 heterologous proteins (ea and eb)
which are involved in substrate degradation
"""

def dydt_mono (t, y, p, substrate):
    ma = y[0]
    rma = y[1]
    ea = y[2] 
    
    mb = y[3] 
    rmb= y[4]
    eb= y[5]
    
    mq= y[6]
    rmq= y[7]
    q= y[8]
    
    mm = y[9]
    rmm = y[10]
    em = y[11]
    
    mt = y[12]
    rmt= y[13]
    et= y[14]
    
    mr = y[15]
    rmr = y[16]
    r = y[17]
    
    si= y[18]
    a= y[19]
    
    mc = y[20]
    rmc = y[21]
    ec = y[22]
    
    z = y[23]
    
    s0 = y[24]
    s1= y[25] 
    s= y[26]
    
    
    #variable constants
    gamma= (p["gmax"]*a/(p["Kgamma"] + a))
    ttrate = ((rmq + rmr + rmt + rmm + rma + rmb + rmc)*gamma)
    
    #growth
    lam = ttrate/p["M"]
    
    #import and metabolism rates
    rate_import = p["v_et"]*et*s/(p["Km_et"] + s)
    rate_metabolism = p["v_em"]*em*si/(p["Km_em"] + si)
    
    #substrate degradation, import & metabolism rates
    
    #if glucose is the primary substrate: 
    if substrate == "glucose":
        ds0dt = 0
        ds1dt = 0      
        dsdt = p["s_in"] - rate_import * (p["Na"] + p["Nb"]) - p["d_s"]*s
    
    
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
            

        ds0dt = p["s0_in"] - rate_ea * (p["Na"] + p["Nb"]) - p["d_s0"]*s0
        ds1dt= rate_ea * (p["Na"] + p["Nb"]) - rate_eb * (p["Na"] + p["Nb"]) - p["d_s1"]*s1
        dsdt = rate_eb * (p["Na"] + p["Nb"]) - rate_import *(p["Na"] + p["Nb"]) - p["d_s"]*s
        
    ####CHANGE FROM EXTRACELLULAR TO SINGLE CELL HERE
        
    #Imported glucose 
    dsidt= rate_import - rate_metabolism - lam*si 
    #a available energy
    dadt= p["ns"]*rate_metabolism - ttrate - lam*a
    
    ###T&T for proteins
    
    #enzyme a
    
    #ma enzyme a mRNA
    dmadt= (p["w_ea"]*a/(p["theta_x"] + a)) + p["ku"]*rma + gamma/p["n_a"]*rma - p["kb"]*r*ma - p["dm"]*ma - lam*ma
    #rma ribo:mRNA complex 
    drmadt= p["kb"]*r*ma - p["ku"]*rma - (gamma/p["n_a"])*rma - lam*rma
    #enzyme alpha (a-amylase)
    deadt= (gamma/p["n_a"])*rma - lam*ea - p["d_ea"] #- p["dg_ea"]*ea add if degradation tag
    
    #enzyme b
    
    #mb enzyme b mRNA
    dmbdt= (p["w_eb"]*a/(p["theta_x"] + a)) + p["ku"]*rmb + gamma/p["n_b"] * rmb - p["kb"]*r*mb - p["dm"]*mb - lam*mb
    #rmb ribo:mRNA complex
    drmbdt= p["kb"]*r*mb - p["ku"]*rmb - (gamma/p["n_b"])*rmb - lam*rmb
    #enzyme b (glucoamylase)
    debdt= (gamma/p["n_b"])*rmb - lam*eb - p["d_ea"] #- p["dg_eb"]*eb add if degradation tag
    
    #housekeeping
    
    #mq housekeeping protein mRNA
    dmqdt= (p["w_q"]*a/(p["theta_x"] + a)) * (1/(1 + (q/p["Kq"])**p["nq"])) + (p["ku"]*rmq) + (rmq*(gamma/p["n_x"])) - (p["kb"]*r*mq) - (p["dm"]*mq) - (lam*mq)
    #rmq ribo:housekeeping mRNA
    drmqdt= (p["kb"]*r*mq) - (p["ku"]*rmq) - (gamma/p["n_x"])*rmq - lam*rmq
    #q housekeeping protein
    dqdt= (gamma/p["n_x"])*rmq - (lam*q)
    
    #metabolic
    
    #mm metabolic enzyme mRNA
    dmmdt= (p["w_m"]*a/(p["theta_x"] + a)) + p["ku"]*rmm + (gamma/p["n_x"])*rmm - p["kb"]*r*mm - p["dm"]*mm - lam*mm
    #rmm Ribo:metabolic enzyme mRNA
    drmmdt= p["kb"]*r*mm - p["ku"]*rmm - (gamma/p["n_x"])*rmm - lam*rmm
    #em enzyme metabolic
    demdt= (gamma/p["n_x"])*rmm - lam*em
    
    #transport
    
    #mt transport enzyme mRNA
    dmtdt= (p["w_t"]*a/(p["theta_x"] + a)) + p["ku"]*rmt + gamma/p["n_x"]*rmt - p["kb"]*r*mt - p["dm"]*mt - lam*mt
    #rmt Ribo:transport enzyme mRNA
    drmtdt= p["kb"]*r*mt - p["ku"]*rmt - (gamma/p["n_x"])*rmt - lam*rmt    
    #et enzyme transport
    detdt= (gamma/p["n_x"])*rmt - lam*et

    #ribosomes

    #mr ribo mRNA
    dmrdt= (p["w_r"]*a/(p["theta_r"] + a)) - (p["kb"]*r*mr) + (p["ku"]*rmr) + rmr*(gamma/p["n_r"]) - (p["dm"]*mr) - (lam*mr)
    #rmr ribo:ribo mRNA
    drmrdt= (p["kb"]*r*mr) - (p["ku"]*rmr) - (gamma/p["n_r"])*rmr - lam*rmr
    #r ribosomes
               #ribsome production     #ribosome unbinding of all mrna complexes                                            #completion of translation for all proteins                                                                                                                                 #binding of free ribosome to mrna                                                                    #dilution
    drdt=  ((gamma/p["n_r"])*rmr) + p["ku"]*rma + p["ku"]*rmb + p["ku"]*rmc + p["ku"]*rmq + p["ku"]*rmm + p["ku"]*rmt + p["ku"]*rmr + rma*(gamma/p["n_a"]) + rmb*(gamma/p["n_b"]) + rmc*(gamma/p["n_c"]) + rmq*(gamma/p["n_x"]) + rmm*(gamma/p["n_x"]) + rmt*(gamma/p["n_x"]) + rmr*(gamma/p["n_r"]) - p["kb"]*r*ma - p["kb"]*r*mb - p["kb"]*r*mc - p["kb"]*r*mq - p["kb"]*r*mm - p["kb"]*r*mt - p["kb"]*r*mr - lam*r

    ####EXTRA PROTEIN - if you want to express a third heterologous protein 

    #enzyme c and associated reaction

    #mx enzyme c mRNA
    dmcdt= (p["w_ec"]*a/(p["theta_x"] + a)) + p["ku"]*rmc + gamma/p["n_x"]*rmc - p["kb"]*r*mc - p["dm"]*mc - lam*mc
    #rmc ribo:mRNA complex 
    drmcdt= p["kb"]*r*mc - p["ku"]*rmc - (gamma/p["n_x"])*rmc - lam*rmc
    #enzyme c (arbitrary enzyme for desirable metabolic pathway)
    decdt= (gamma/p["n_c"])*rmc - lam*ec - p["dg_ec"]*ec    

    #reaction converting "energy" (ATP) into a product
    
    dzdt= (p["v_ec"] * ec * a / (p["Km_ec"] + a)) - (lam * z)

    return (dmadt, drmadt, deadt, dmbdt, drmbdt, debdt, dmqdt, drmqdt, dqdt, dmmdt, drmmdt, demdt, dmtdt, drmtdt, detdt, dmrdt, drmrdt, drdt, dsidt, dadt, dmcdt, drmcdt, decdt, dzdt, ds0dt, ds1dt, dsdt)


