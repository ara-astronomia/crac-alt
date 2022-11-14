import configparser as conf
import PySimpleGUI as sg
import os.path #vediamo se sarà necessario
import numpy as np
import matplotlib.pyplot as plt
import datetime
import astropy.units as u
from astropy.visualization import astropy_mpl_style, quantity_support
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astroplan.plots import plot_altitude
from astroplan import Observer
from astropy.coordinates import get_sun
from astropy.coordinates import get_moon
from pytz import timezone
plt.style.use(astropy_mpl_style)
quantity_support()

config = conf.ConfigParser()
config.read('config.ini')

###read value from config.ini
lat = config.get("geografica", "lat")
lon = config.get("geografica", "lon")
elev = config.get("geografica", "elev")
utc_offset = config.get("observer", "utcoffset")
testo_target= config.get("tooltip", "testo_target")
testo_ar= config.get("tooltip", "testo_ar")
testo_dec= config.get("tooltip", "testo_dec")
ref_midnight=config.get("observer", "reference_midnight")
description=config.get("observer", "description")

today=datetime.date.today()

#calcola i parmaetri per la geografia
frasso_sabino = EarthLocation.from_geodetic(float(lon)*u.deg, float(lat)*u.deg, float(elev)*u.m, ellipsoid=None)
utcoffset = +((int(utc_offset))*u.hour)  # UTC Time?

#converte le coordinate inserite e calcaola il valore Altitude per il target e plotta i risultati
def skycoord(ar, dec, time_obs, name_target, local_time):
    ar = ar
    dec = dec
    
    time_obs= Time(time_obs +' '+ ref_midnight, scale='utc')
    print (time_obs)

    if local_time is True:
        time_obs=time_obs+utcoffset
        print (time_obs)
      
    target =SkyCoord((ar +' ' +dec), frame=FK5, unit=(u.hourangle, u.deg), obstime=time_obs)
    targetaltaz = target.transform_to(AltAz(obstime=time_obs,location=frasso_sabino))
    print(f"TARGET's Altitude = {targetaltaz.alt:.2}")
    coord=(ar +' '+dec)
  
    #calcola la mezzanotte
    midnight = time_obs #- utcoffset
    delta_midnight = np.linspace(-8, 24, 100)*u.hour

    #calcola il deltamidnight per l'alteza del sole
    delta_midnight = np.linspace(-8, 24, 1000)*u.hour
    times_today_up_tomorrow = midnight + delta_midnight
    
    frame_today_up_tomorrow = AltAz(obstime=times_today_up_tomorrow, location=frasso_sabino)
    sunaltazs_today_up_tomorrow = get_sun(times_today_up_tomorrow).transform_to(frame_today_up_tomorrow)

    #calcola la posizione della luna
    moon_today_up_tomorrow = get_moon(times_today_up_tomorrow, location=frasso_sabino)
    moonaltazs_today_up_tomorrow = moon_today_up_tomorrow.transform_to(frame_today_up_tomorrow)

    # Find the alt,az coordinates of target at those same times:
    targetaltazs_today_up_tomorrow = target.transform_to(frame_today_up_tomorrow)


    plt.rc('font', size=8)
    plt.plot(delta_midnight, sunaltazs_today_up_tomorrow.alt, color='r', label='Sun')
    plt.plot(delta_midnight, moonaltazs_today_up_tomorrow.alt, color=[0.75]*4, ls='--', label='Moon')
    plt.scatter(delta_midnight, targetaltazs_today_up_tomorrow.alt,
                c=targetaltazs_today_up_tomorrow.az, label=name_target+' '+coord, lw=0, s=8,
                cmap='viridis')
    plt.fill_between(delta_midnight, 0*u.deg, 90*u.deg,
                    sunaltazs_today_up_tomorrow.alt < -0*u.deg, color='0.8', zorder=0)
    plt.fill_between(delta_midnight, 0*u.deg, 90*u.deg,
                    sunaltazs_today_up_tomorrow.alt < -18*u.deg, color='k', zorder=0)
    #plt.colorbar().set_label('Azimuth [deg]')

    plt.legend(loc='upper left')
    plt.xlim(-8*u.hour, 8*u.hour)
    plt.xticks((np.arange(18)*1-9)*u.hour)
    plt.ylim(0*u.deg, 90*u.deg)
    plt.yticks((np.arange(18)*5)*u.deg)
    plt.xlabel('Hours from UTC Midnight @ '+str(time_obs) +'\n '+' in '+ str(description)+' lat '+str(lat) +' '+'lon '+ str(lon))
    plt.ylabel('Altitude [deg]')
    plt.show()

# Create the Window
sg.theme('Dark')
layout = [
            [sg.Text('Observatory', size=10), sg.Text(description, size=25)],
            [sg.Text('Lat.', size=10), sg.Text(lat , size=18)],
            [sg.Text('Lon.', size=10), sg.Text(lon , size=18)],
            [sg.Text('UTC Offset', size=10), sg.Text(utcoffset , size=8)],
            [sg.Checkbox('ora locale', key='_LOCAL_TIME_')],
            [sg.Text('Time:', size=8), sg.Input(today, key='iel', size=15), sg.CalendarButton('Cambia la data', target='iel',format='%Y-%d-%m ')],
            [sg.Text("target", size=8), sg.InputText(key='_TARGET_', size=15, tooltip=testo_target)],
            [sg.Text("AR ", size=8), sg.InputText(key='_AR_KEY_', size=15, tooltip=testo_ar)],
            [sg.Text("Dec", size=8), sg.InputText(key='_DEC_KEY_', size=15, tooltip=testo_dec)],        
            [sg.Button('Calcola'), sg.Button('Exit')],
            
          ]


window = sg.Window('CRaCAlt by ARA', layout)
# Create the event loop
while True:
    event, values = window.read()
    if event in ('Calcola'):
        print('ora calcolo')
        ar=values['_AR_KEY_']
        dec=values['_DEC_KEY_']
        name_target=values['_TARGET_']
        time_obs=values['iel']
        local_time=values['_LOCAL_TIME_']
        print(time_obs, local_time)
                      
        if not name_target:
            sg.popup('il campo target è vuoto, verrà impostato il nome di default Target')
            name_target='Target'
            
        skycoord(ar, dec, time_obs, name_target, local_time)
       
    if event in (None, 'Exit'):
        # User closed the Window or hit the Cancel button
        break
    print(f'Event: {event}')
    #print(str(values))

window.close()