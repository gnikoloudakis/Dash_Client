'''
Author: asideris
Email:	 asiderisa@gmail.com
Created: Fri Feb 28 13:01:21 2014
'''

import sys
import m3u8
import time
import json
import math
import thread
import requests
from lxml import etree
from collections import OrderedDict


class ManifestParser(object):
    """This class defines the
    manifests' parsing operations"""

    def __init__(self):
        "initialise's the class"
        # print inspect.getfile(inspect.currentframe())
        # print "Object " + self.__class__.__name__ + " initialised"
        pass

    def manifest_type(self, manifest):
        """
        INPUT: A manifest's uri (:string)
        OUTPUT: Manifest's type (:string)
        USES: -
        DESC: This method discovers and returns the manifest type (i.e HLS,
        MPEG-DASH)
        """

        mtype = ""
        try:
            if manifest.endswith('.mpd'):
                mtype = 'MPEG-DASH'
            if manifest.endswith('.m3u8'):
                mtype = 'HLS'
            return mtype
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None

    def manifest_variancy(self, manifest):
        """
        INPUT: A manifest's uri (:string)
        OUTPUT: Manifest's variancy (:boolean)
        USES: -
        DESC: This method finds out if a manifest has more than one
        quality profiles to offer.
        """

        variancy = False
        try:
            mType = self.manifest_type(manifest)

            if mType == "MPEG-DASH":
                parser = etree.XMLParser(ns_clean=False,
                                         remove_blank_text=True)
                xmlDoc = etree.parse(manifest, parser)
                root = xmlDoc.getroot()

                nsUri = root.nsmap.items()[0][1]

                # set Xpath search string for representations
                ssRep = etree.XPath("//p:Representation",
                                    namespaces={'p': nsUri})

                # Find all mpeg-dash representations
                repList = (ssRep(root))

                if len(repList) > 1:
                    variancy = True
                    # print len(repList)
                    pass
                pass

            if mType == "HLS":
                mf = m3u8.load(manifest)

                # checks if it points to other m3u8 files (is variant)
                if mf.is_variant:
                    variancy = True
                    pass
                pass

            return variancy
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None

    def hls_parser(self, manifest):
        """
        INPUT: A manifest's uri (:string)
        OUTPUT: Manifest's playlist (:OrderedDict)
        USES: hls_Variant_parser(), hls_nonVariant_parser()
        DESC: This method initiates the parsing of a HLS
        manifest and returns a playlist for the DASH
        Client emulator. The playlist is stored in an
        OrderedDictionary having the following format:
        {qualityProfile:{manifest:'',rate:'',resolution:''\
        ,codecs:'',segments:{id:{segment:'',duration:''}}}
        """

        playlist = OrderedDict()
        try:

            mf = m3u8.load(manifest)

            # checks if it points to other m3u8 files (is variant)
            if mf.is_variant:
                return self.hls_variant_parser(mf)
            else:
                playlist[1] = OrderedDict(
                    {'manifest': manifest, 'rate': 'NA',
                     'resolution': 'NA', 'codecs': 'NA',
                     'segments': self.hls_nonVariant_parser(mf)})
                # print playlist
                return playlist
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None
        finally:
            print
            'Exiting hls_parser'

    def hls_variant_parser(self, variantManifest):
        """
        INPUT: A variant manifest (:m3u8.model.M3U8)
        OUTPUT: Manifest's playlist (:OrderedDict)
        USES: hls_nonVariant_parser()
        DESC: This method parses a variant HLS manifest
        and returns a playlist for the DASH Client emulator.
        The playlist is stored in an OrderedDictionary having
        the following format:
        {qualityProfile:{manifest:'',rate:'',resolution:''\
        ,codecs:'',segments:{id:{segment:'',duration:''}}}
        """

        playlist = OrderedDict()
        try:
            index = 1
            for profile in variantManifest.playlists:
                playlist[index] = OrderedDict(
                    {'manifest': profile.uri,
                     'rate': profile.stream_info.bandwidth,
                     'resolution': profile.stream_info.resolution,
                     'codecs': profile.stream_info.codecs,
                     'segments': self.hls_nonVariant_parser(
                         m3u8.load(profile.uri))})

                index = index + 1

            return playlist
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None

    def hls_nonVariant_parser(self, nonVariantManifest):
        """
        INPUT: A non variant manifest (:m3u8.model.M3U8)
        OUTPUT: Manifest's playlist (:OrderedDict)
        USES: -
        DESC: This method parses a non variant HLS manifest
        and returns the segments' location and duration.
        The playlist is stored in an
        OrderedDictionary having the following format:
        {id:{segment:'',duration:''}}
        """

        try:
            index = 1
            segments = OrderedDict()

            for segment in (nonVariantManifest).segments:
                location = str(segment).split()[1]
                duration = str(segment).split()[0].split(':')[1].split(',')[0]
                segments[index] = {'segment': location, 'duration': duration}
                index = index + 1
            return segments
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None

    def mpegDash_parser(self, manifest):
        """
        INPUT: A manifest's uri (:string)
        OUTPUT: Manifest's playlist (:OrderedDict)
        USES: mpegDash_representations_parser()
        DESC: This method initiates the parsing of a
        MPEG-DASH manifest and returns a playlist for
        the DASH Client emulator. The playlist is stored
        in an OrderedDictionary having the following format:
        {representation:{baseUrl:'',rate:'',resolution:''\
        ,codecs:'',segments:{id:{segment:'',duration:''}}}
        """

        try:
            parser = etree.XMLParser(ns_clean=False, remove_blank_text=True)
            xmlDoc = etree.parse(manifest, parser)
            root = xmlDoc.getroot()

            # get namespace
            nsUri = root.nsmap.items()[0][1]
            # set Xpath search string for Base URL
            ssBaseURL = etree.XPath("//p:BaseURL", namespaces={'p': nsUri})
            # set Xpath search string for representations
            ssRep = etree.XPath("//p:Representation", namespaces={'p': nsUri})
            # Find Base URL
            baseURL = (ssBaseURL(root))[0].text
            # Find all mpeg-dash representations
            repList = (ssRep(root))

            if len(repList) == 0:
                print("No representations found")
                return None
            else:
                playlist = OrderedDict()
                for rep in repList:
                    playlist[int(rep.attrib['id'])] = OrderedDict(
                        {'baseUrl': baseURL,
                         'rate': rep.attrib['bandwidth'],
                         'resolution': rep.attrib['width']
                                       + 'x' + rep.attrib['height'],
                         'codecs': rep.attrib['codecs'],
                         'segments': self.mpegDash_representations_parser(
                             rep, baseURL)})

                    # print type(rep), type(baseURL)

                return playlist

        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None
        finally:
            print('Exiting mpegDash_parser')

    def mpegDash_representations_parser(self, representation, baseURL):
        """
        INPUT: A mpeg-DASH representation, a baseURL
        (:lxml.etree._Element, :string)
        OUTPUT: Manifest's playlist (:OrderedDict)
        USES: -
        DESC: This method parses a MPEG-DASH representation
        and returns the segments' location and duration.
        The playlist is stored in an
        OrderedDictionary having the following format:
        {id:{segment:'',duration:''}}
        """

        try:
            index = 1
            segments = OrderedDict()

            # get namespace
            nsUri = representation.nsmap.items()[0][1]

            # set Xpath search string for Segments list
            fSgmList = etree.XPath(".//p:SegmentList",
                                   namespaces={'p': nsUri})
            # Find Segments list
            sgmList = (fSgmList(representation))
            # set Xpath search string for Segments URL list
            fSgmURLList = etree.XPath(".//p:SegmentURL",
                                      namespaces={'p': nsUri})
            # # Find Segments URL list
            sgmURLList = (fSgmURLList(representation))

            # get Segment  duration
            duration = float(
                sgmList[0].attrib['duration']
            ) / float(sgmList[0].attrib['timescale'])

            for segment in sgmURLList:
                segments[index] = {'segment': baseURL + segment.attrib['media'],
                                   'duration': duration}
                index = index + 1
            return segments
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None


class QualityProfiles(object):
    """
    This class provides and coordinates
    all the quality profiles
    selection related operations
    """

    def __init__(self):
        "Initialises class"
        # print "Object " + self.__class__.__name__ + " initialised"
        pass

    def min_qualityProfile(self, plist):
        """
        INPUT: a manifest's playlist  (:OrderedDict)
        OUTPUT: min quality profile's id (:int)
        USES: -
        DESC: This method finds and returns minimum quality profile's id.
        """
        try:
            minProfile = None
            ## set min value to the first playlist's rate ##
            minRate = float(plist.items()[0][1]['rate'])
            # print minRate
            for k, v in plist.items():
                if float(v['rate']) <= minRate:
                    minRate = float(v['rate'])
                    minProfile = k
            return minProfile
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None

    def get_profile_segment_count(self, profile, plist):
        """
        INPUT: a profile's id (:int), a manifest's playlist  (:OrderedDict)
        OUTPUT: min quality profile's id (:int)
        USES: -
        DESC: This method finds and returns minimum quality profile's id.
        """
        try:
            return len(plist.get(profile).get('segments'))
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None

    def get_profile_rate(self, profile, plist):
        """
        INPUT: a profile's id (:int),
        manifest's playlist  (:OrderedDict)
        OUTPUT: the profile's rate (:float)
        USES: -
        DESC: This method returns the profile's rate.
        """
        try:
            return float(plist.get(profile).get('rate'))
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None

    def get_profile_id(self, rate, plist):
        """
        INPUT: a profile's rate (:float),
        manifest's playlist  (:OrderedDict)
        OUTPUT: the profile's id (:int)
        USES: -
        DESC: This method returns a profile's id.
        """
        try:
            for k, v in plist.items():
                if float(v['rate']) == rate:
                    return k
            print
            "No matching rate"
            return -1
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None

    def get_profiles_id(self, plist):
        """
        INPUT: a manifest's playlist  (:OrderedDict)
        OUTPUT: a list with the profiles id (:list)
        USES: -
        DESC: This method returns the playlist's profiles id.
        """
        try:
            idList = []
            for k, v in plist.items():
                idList.append(k)
            return idList
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None


class Logic(object):
    """
    This class provides the
    available HAS logics
    """

    def __init__(self):
        "Class initialisation"
        # print "Object " + self.__class__.__name__+" initialised"

    def first_time_l1(self, plist):
        """
        INPUT: a manifest's playlist  (:OrderedDict)
        OUTPUT: min quality profile's id (:int)
        USES: min_qualityProfile()
        DESC: This logic is used select the minimum
        quality profile of a playlist.
        """
        try:
            return QualityProfiles().min_qualityProfile(plist)
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None

    def logic_1(self, sRate, profile, plist):
        """
        INPUT: the segment's download rate (:float),
        current profile's id (:int), a manifest's
        playlist  (:OrderedDict)
        OUTPUT: switch indicator (:boolean), select indicator (:boolean) and
        selected profile's id
        USES: -
        DESC: This method compares a segment's download rate (d1)
        with the active profile's encoding rate (p1) and if d1 != p1
        selects the max profile having rate p2 <= d1. 
        If no profile is available it returns the minimum rate profile 
        """
        try:
            switchP = False
            pRates = []
            newProfile = profile

            pRate = float(plist.get(profile).get('rate'))
            # print pRate
            if pRate != sRate:
                switchP = True

            if switchP:
                ''' get playlist's profiles id '''
                idList = QualityProfiles().get_profiles_id(plist)

                for id in idList:
                    pRates.append(float(plist.get(id).get('rate')))

                ''' get profiles fitting to our rate threshold '''
                avlRates = [r for r in pRates if r <= sRate]

                ''' Check and find are available profile rates '''
                if avlRates:
                    for id in idList:
                        if float(plist.get(id).get('rate')) == max(avlRates):
                            if id != profile:
                                newProfile = id
                            break
                else:
                    newProfile = QualityProfiles().min_qualityProfile(plist)

            return [switchP, newProfile != profile, newProfile]
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None

    def logic_2(self, capacity, estRate, profile, plist):
        """
        INPUT: the bottleneck link's capacity (:float), the estimated bit rate (:float),
        current profile's id (:int), a manifest's
        playlist  (:OrderedDict)
        OUTPUT: switch indicator (:boolean), select indicator (:boolean) and
        selected profile's id
        USES: -
        DESC: This method provides paper's
        'Improving Fairness in QoS and QoE domains for Adaptive Video Streaming'
        switching/selecting logic.
        """
        try:
            switchP = False
            pRates = []
            newProfile = profile
            congThr = 0.95 * capacity
            # burstThr = 1.1 * capacity

            switchP = estRate < congThr or estRate > capacity

            if switchP:
                ''' get currrent's profile rate '''
                pRate = float(plist.get(profile).get('rate'))

                ''' get playlist's profiles id '''
                idList = QualityProfiles().get_profiles_id(plist)

                for id in idList:
                    pRates.append(float(plist.get(id).get('rate')))

                if estRate < congThr:
                    print
                    'Xa'
                    ''' get profiles fitting to our rate threshold '''
                    avlRates = [r for r in pRates if r <= estRate]
                    ''' Check and find the available profile rates '''
                    if avlRates:
                        suppRate = max(avlRates)
                        newProfile = QualityProfiles().get_profile_id(suppRate, plist)
                    else:
                        newProfile = QualityProfiles().min_qualityProfile(plist)

                if estRate > capacity:
                    print
                    'Xo'
                    ''' get profiles fitting to our rate threshold '''
                    avlRates = [r for r in pRates if r < pRate]
                    ''' Check and find are available profile rates '''
                    if avlRates:
                        suppRate = max(avlRates)
                        newProfile = QualityProfiles().get_profile_id(suppRate, plist)
                    else:
                        newProfile = QualityProfiles().min_qualityProfile(plist)

            return [switchP, newProfile != profile, newProfile]
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None

    def logic_3(self, dwlTime, sgmDur, profile, plist):
        """
        INPUT: the segment's download time (:float), the segment's duration (:float),
        current profile's id (:int), a manifest's
        playlist  (:OrderedDict)
        OUTPUT: switch indicator (:boolean), select indicator (:boolean) and
        selected profile's id
        USES: -
        DESC: This method provides Adobe's OSMF logic as paper
        'QDASH: A QoE-aware DASH system' indicates.
        """
        try:
            switchP = False
            pRates = []
            newProfile = profile

            rDwl = sgmDur / float(dwlTime)

            '''get minimum profile id'''
            minP = QualityProfiles().min_qualityProfile(plist)
            ''' get minimum profile rate '''
            minPrate = float(plist.get(minP).get('rate'))
            ''' get currrent's profile rate '''
            pRate = float(plist.get(profile).get('rate'))

            ''' get playlist's profiles id '''
            idList = QualityProfiles().get_profiles_id(plist)
            for id in idList:
                pRates.append(float(plist.get(id).get('rate')))

            ''' get previous profile rate '''
            if pRate == minPrate:
                prvRate = minPrate
            else:
                prvRate = max([r for r in pRates if r < pRate])

            # print '{0} {1} {2}'.format(rDwl, sgmDur, dwlTime)
            if rDwl < 1:

                if pRate > minPrate:
                    switchP = True
                    if rDwl < (prvRate / float(pRate)):
                        newProfile = minP
                    else:
                        newProfile = QualityProfiles().get_profile_id(prvRate, plist)
            else:
                if pRate < max(pRates):
                    if rDwl >= (prvRate / float(pRate)):
                        switchP = True
                        currRate = pRate
                        while True:
                            nxtRate = min([r for r in pRates if r > pRate])
                            if nxtRate == max(pRates):
                                nxtNxtRate = nxtRate
                            else:
                                nxtNxtRate = min([r for r in pRates if r > nxtRate])
                            if nxtRate == max(pRates) or rDwl < (nxtNxtRate / currRate):
                                newProfile = QualityProfiles().get_profile_id(nxtRate, plist)
                                break
                            pRate = nxtRate

            return [switchP, newProfile != profile, newProfile]
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None

    def logic_4(self, sgmRate, sgmDur, pBuffer, profile, plist):
        """
        INPUT: the segment's download rate (:float), the segment's duration (:float),
        playout buffer current size (:float), current profile's id (:int), a manifest's
        playlist  (:OrderedDict)
        OUTPUT: switch indicator (:boolean), select indicator (:boolean) and
        selected profile's id
        USES: -
        DESC: This method provides the logic of
        'QDASH: A QoE-aware DASH system' paper.
        """
        try:
            switchP = False
            pRates = []
            newProfile = profile

            '''get minimum profile id'''
            minP = QualityProfiles().min_qualityProfile(plist)
            ''' get minimum profile rate '''
            minPrate = float(plist.get(minP).get('rate'))
            ''' get currrent's profile rate '''
            pRate = float(plist.get(profile).get('rate'))

            ''' get playlist's profiles id '''
            idList = QualityProfiles().get_profiles_id(plist)
            for id in idList:
                pRates.append(float(plist.get(id).get('rate')))

            ''' get supported profile's rate '''
            if minPrate > sgmRate:
                suppRate = minPrate
            else:
                suppRate = max([r for r in pRates if r <= sgmRate])

            if suppRate < pRate:
                switchP = True
                suppId = QualityProfiles().get_profile_id(suppRate, plist)
                currId = QualityProfiles().get_profile_id(pRate, plist)
                if math.fabs(currId - suppId) > 1:
                    sFrag = suppRate * sgmDur
                    tBuff = pBuffer / (1 - (suppRate / sFrag))
                    nFrag = tBuff * (suppRate / sFrag)
                    if nFrag > 0:
                        if suppRate == max(pRates):
                            nxtSuppRate = suppRate
                        else:
                            nxtSuppRate = min([r for r in pRates if r > suppRate])
                        newProfile = QualityProfiles().get_profile_id(nxtSuppRate, plist)
                    else:
                        newProfile = QualityProfiles().get_profile_id(suppRate, plist)
                else:
                    newProfile = QualityProfiles().get_profile_id(suppRate, plist)
            else:
                if suppRate != pRate:
                    switchP = True
                newProfile = QualityProfiles().get_profile_id(suppRate, plist)

            return [switchP, newProfile != profile, newProfile]
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None

    def logic_5(self, sgmRate, sgmDur, pBuffer, profile, plist):
        """
        INPUT: the segment's download rate (:float), the segment's duration (:float),
        playout buffer current size (:float), current profile's id (:int), a manifest's
        playlist  (:OrderedDict)
        OUTPUT: switch indicator (:boolean), select indicator (:boolean) and
        selected profile's id
        USES: -
        DESC: This method provides the logic of
        'Adaptation Algorithm for Adaptive Streaming over HTTP' paper.
        """
        try:
            switchP = False
            pRates = []
            newProfile = profile

            '''get minimum profile id'''
            minP = QualityProfiles().min_qualityProfile(plist)
            ''' get minimum profile rate '''
            minPrate = float(plist.get(minP).get('rate'))
            ''' get currrent's profile rate '''
            pRate = float(plist.get(profile).get('rate'))

            ''' get playlist's profiles id '''
            idList = QualityProfiles().get_profiles_id(plist)
            for id in idList:
                pRates.append(float(plist.get(id).get('rate')))

            return [switchP, newProfile != profile, newProfile]
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None

    def est_bw_l2(self, sRate, avgRate, varRate):
        """
        INPUT: the segment's download rate (:float),
        the avg Rate (:float), the bit rate variancy (:float)
        OUTPUT: the estimated bit rate (:float), the next avg Rate (:float),
        the next bit rate variancy (:float)
        USES: -
        DESC: This method evaluates the achievable bit rate
        as the 'DISTRIBUTED & ADAPTIVE HTTP STREAMING' paper indicates.
        """
        try:
            a = 1.0 / 16
            b = 1.0 / 8
            c = 1.0

            nextAvgRate = (1 - a) * avgRate + (a * sRate)
            diff = math.fabs(sRate - avgRate)
            nextVarRate = (1 - b) * varRate + (b * diff)
            estRate = avgRate - (c * varRate)

            return [estRate, nextAvgRate, nextVarRate]
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None

    def est_bw_l3(self, sRate, prevRate):
        """
        INPUT: the segment's download rate (:float),
        the previous Rate (:float)
        OUTPUT: the estimated bit rate (:float)
        USES: -
        DESC: This method evaluates the achievable bit rate
        as the 'Improving Fairness in QoS and QoE domains for Adaptive Video Streaming' paper indicates.
        """
        try:
            d = 0.8

            return (d * prevRate) + ((1 - d) * sRate)
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None


class Segments(object):
    """
    This class provides segment
    related  operations
    """

    def __init__(self):
        "Class initialisation"
        # print "Object " + self.__class__.__name__+" initialised"

    def get_segment(self, uri):
        """
        INPUT: A segment's uri (:string)
        OUTPUT: download statistics (:OrderedDict)
        USES: -
        DESC: This method returns the download statistics
        of each segment in a dictionary, which has the
        following format:
        {dwlBytes:'Bytes',dwlTime:'Seconds',dwlRate:'bps',
        'startDwlTime:'Seconds',endDwlTime':'Seconds'}
        """

        try:
            statsDic = OrderedDict()
            startTime = round(time.time(), 6)
            r = requests.get(uri)
            endTime = round(time.time(), 6)

            ''' Elapsed Time in seconds '''
            elapsedTime = endTime - startTime
            # print elapsedTime
            ''' Converts Bps to bps '''
            dwlRate = ((len(r.content) / elapsedTime) * 8)
            # print ("finish downloading at "+str(dwlRate)+" bps")

            statsDic['startDwlTime'] = startTime
            statsDic['endDwlTime'] = endTime
            statsDic['dwlBytes'] = len(r.content)  # Unit=Bytes
            statsDic['dwlTime'] = elapsedTime  # Unit=Seconds
            statsDic['dwlRate'] = dwlRate  # Unit=bps
            return statsDic
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None

    def get_segment_uri(self, index, profile, playlist):
        """
        INPUT: A segment's index (:int), a profile's index (:int),
        a manifest's playlist (:OrderedDict)
        OUTPUT: A segment's uri (:string)
        USES: -
        DESC: This method returns the uri of a segment
        """
        try:

            uri = playlist.get(
                profile).get('segments').get(index).get('segment')
            return uri
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            print
            "Check the playlist's keys"
            return None

    def get_segment_duration(self, index, profile, playlist):
        """
        INPUT: A segment's index (:int), a profile's index (:int),
        a manifest's playlist (:OrderedDict)
        OUTPUT: A segment's duration (:float)
        USES: -
        DESC: This method returns the duration of a segment
        """
        try:

            dur = playlist.get(
                profile).get('segments').get(index).get('duration')
            return float(dur)
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            print
            "Check the playlist's keys"
            return None


class Statistics(object):
    """
    This class provides the
    operations related with
    Statistics Info
    """

    def __init__(self):
        "Initialises the class"

        # print "Object " + self.__class__.__name__ + " initialised"

    def write_stats(self, statsDic, statsFile):
        """
        INPUT: segment Stats (:OrderedDict),
        a stats File (:string)
        OUTPUT: -
        USES: -
        DESC: This method saves/updates the stats file
        with the info contained in the segment stats dic.
        """
        try:
            # with open(statsFile) as f:
            #     stats = json.load(f)

            # stats.update(statsDic)
            with open(statsFile, 'w') as f:
                json.dump(statsDic, f, indent=4, sort_keys=True)

            pass
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))

    def prepare_segment_stats(self, pRate, pid, sgmUri, sgmDur):
        """
        INPUT: profile's rate (:float), profile's id (:int),
        segment's uri (:string) and segment's duration (:float)
        OUTPUT: stats regarding segments (:dict)
        USES: -
        DESC: This method  returns the stats regarding
        segments.
        """
        try:
            sgmInfo = {'profile_rate': pRate,
                       "profile_id": pid,
                       "segment_uri": sgmUri,
                       "segment_duration": sgmDur}
            return sgmInfo
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))

    def prepare_logic_stats(self, swFlag, scFlag, estBw):
        """
        INPUT: Switch indicator flag (:boolean),
        select indication flag (:boolean), estimated BW (:float)
        OUTPUT: stats regarding logic (:dict)
        USES: -
        DESC: This method returns the stats regarding the
        switching/selection logic.
        """
        try:
            logicInfo = {"Switch": swFlag,
                         "Select": scFlag,
                         "Estimated_BW": estBw}
            return logicInfo
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))

    def prepare_pBuffer_stats(self, minPBuf, maxPBuf, currPBuf):
        """
        INPUT: min playout buffer size (:int),
        max playout buffer size (:int), current buffer
        size (:OrderedDict), playout buffer stall (:float)
        OUTPUT: stats regarding playout buffer (:dict)
        USES: -
        DESC: This method returns the stats regarding the
        playout buffer.
        """
        try:
            # pBufInfo = OrderedDict({"min_playout_buffer": minPBuf,
            #                         "max_playout_buffer": maxPBuf,
            #                         "curr_playout_buffer": currPBuf})
            pBufInfo = OrderedDict({"min_playout_buffer": minPBuf,
                                    "max_playout_buffer": maxPBuf,
                                    "curr_playout_buffer": {
                                        "pBuf_1": currPBuf['pBuf_1'],
                                        "pBuf_2": currPBuf['pBuf_2'],
                                        "pBuf_3": currPBuf['pBuf_3']}})
            return pBufInfo
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))

    def prepare_all_stats(self, statsDic, dwlStats, idleDwlTime,
                          sgmInfo, logicInfo, pBufInfo, sgmIndex):
        """
        INPUT: statsDic (:OrderedDict), dwl stats (:dict),
        idle download time (:float), segment's info (:Dict),
        logic's info (:Dict), playout buffer's info (:Dict)
        and segment's index (:int)
        OUTPUT: updated stats  (:dict)
        USES: -
        DESC: This method  returns all the stats in one dictionary.
        """
        try:

            statsDic[sgmIndex] = {
                "segment_info": sgmInfo,
                "logic_info": logicInfo,
                "dwlStats": dwlStats,
                "idleDwlTime": idleDwlTime,
                "playout_buffer_info": pBufInfo}

            # print "\nstats: %s" % pBufInfo
            # statsDic[sgmIndex] = ({"playout_buffer_info": pBufInfo})

            # print statsDic
            return statsDic
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))


class PlayoutBuffer(object):
    """
    This class provides the
    operations related with
    the playout Buffer operations
    """

    def __init__(self):
        "Initialises the class"

        # print "Object " + self.__class__.__name__ + " initialised"

    def check_max_buffer_for_init_buf_state(self, pBuf, maxPBuf, sgmDur):
        """
        INPUT: the playout buffer size (:float), max playout buffer size (:int)
        and segment duration (:float)
        OUTPUT: -
        USES: -
        DESC: This method checks for playout buffer overflow when in
        initial buffering state
        """
        try:
            if (pBuf + sgmDur) > maxPBuf:
                sys.exit("Max Playout Buffer overflow. Please change " +
                         "max threshold or reduce initial " +
                         "buffering segments")
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))

    def minThreshold_buffer_l1(self, pBuf, maxPBuf, minPBuf, sgmDur):
        """
        INPUT: the playout buffer size (:float), max playout
        buffer size (:int), min playout buffer size (:int)
        and segment duration (:float)
        OUTPUT: -
        USES: -
        DESC: This method pauses the dwonlaod operations until the
        playout buffer becomes <= min playout buffer threshold
        """
        try:
            startTime = round(time.time(), 6)

            ''' Check If min Playout Buffer + duration exceeds
            max playout buffersize. Fixing by changing
            temporary the min playout buf size to the value
            new val = maxPbuf - selected segment's duration '''
            if (maxPBuf - minPBuf) < sgmDur:
                minThr = minPBuf - sgmDur
                if minThr < 0:
                    minThr = 0
            else:
                minThr = minPBuf

            while True:
                elapsedTime = round(time.time(), 6) - startTime
                if (pBuf - elapsedTime) <= minThr:
                    return elapsedTime
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))

    def maxThreshold_buffer_l1(self, pBuf, maxPBuf, sgmDur, elapsedTime):
        """
        INPUT: the playout buffer size (:float), max playout
        buffer size (:int), segment duration (:float) and
        the waiting time for reaching the min playout
        buffer threshold  (:float)
        OUTPUT: -
        USES: -
        DESC: This method checks for playout buffer overflow when in
        playout state
        """
        try:
            increase = sgmDur - elapsedTime

            return True if (pBuf['pBuf_3'] + increase) <= maxPBuf else False
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None

    def calc_curr_playout_buffer_for_init_buf_state(self, pBuf, sgmDur):
        """
        INPUT: the current playout buffer (:Dict),
        and segment duration (:float)
        OUTPUT: the updated current playout buffer (:Dict)
        USES: -
        DESC: This method updates the playout buffer for the initial
        buffering state
        """
        try:
            ''' calculate playout buffer size for before segment download '''
            pBuf['pBuf_1'] = pBuf['pBuf_3']
            ''' calculate playout buffer size for after segment download
            not including the downloaded segment '''
            pBuf['pBuf_2'] = pBuf['pBuf_3']
            ''' calculate playout buffer size for after segment download
            including the downloaded segment '''
            pBuf['pBuf_3'] += sgmDur
            return pBuf
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))

    def calc_curr_playout_buffer(self, pBuf, sgmDur, idleDwlTime, sgmDwlTime):
        """
        INPUT: the current playout buffer (:Dict),
        segment duration (:float), idle download time (:float)
        and segment download time (:float)
        OUTPUT: the updated current playout buffer (:Dict)
        USES: -
        DESC: This method updates the playout buffer
        """
        try:
            ''' calculate playout buffer size
            for before segment download '''
            pBuf['pBuf_1'] = pBuf['pBuf_3'] - idleDwlTime

            ''' fill playout buffer size after segment download
            not including the downloaded segment '''
            pBuf['pBuf_2'] = pBuf['pBuf_1'] - sgmDwlTime

            '''fill playout buffer size after segment download
            including the downloaded segment '''
            if pBuf['pBuf_2'] >= 0:
                pBuf['pBuf_3'] = pBuf['pBuf_2'] + sgmDur
            else:
                pBuf['pBuf_3'] = sgmDur

            return pBuf
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))


class DashClientEmulator(object):
    """
    This is the root class.
    It will provide and coordinate
    all the DASH related operations
    """

    def __init__(self):
        "Initialises root class"
        print("Object " + self.__class__.__name__ + " initialised")

    def parse_operation(self):
        """
        INPUT: self
        OUTPUT: a tuple of the following form
        [manifest's playlist (:OrderedDict), variancy (:boolean)]
        USES: ManifestParser class methods
        DESC: This method initiates the parse operations and returns a
        a playlist and a variancy indicator
        """
        try:
            ''' get manifest's uri '''
            mUri = self.manifest

            plist = None
            mParser = ManifestParser()

            mType = mParser.manifest_type(mUri)
            # print mType
            if mType == 'MPEG-DASH':
                print(mType)
                plist = mParser.mpegDash_parser(mUri)
                pass

            if mType == "HLS":
                print(mType)
                plist = mParser.hls_parser(mUri)
                pass

            pass
            return [plist, mParser.manifest_variancy(mUri)]
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None

    def info(self):
        """
        INPUT: self
        OUTPUT:-
        USES: -
        DESC: this method print some statistics per downloaded segment.
        """
        try:
            print("\nSegment %s" % self.sgmIndex)

            print("Segment's uri: %s" % self.sgmUri)

            print("Current  Dwl Rate {0} bps --- Estimated Rate {1} bps".format(self.sgmRate, self.estRate))

            print("current Profile: {0} -- next Profile: {1}".format(self.previousProfile, self.currentProfile))

            print("current playout buffer is %s" % self.curPlayBuf)

            print("Idle download time: %s sec" % self.idleDwlTime)

        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None

    def write_stats_wrapper(self):
        """
        INPUT: self
        OUTPUT:-
        USES: Statistics class methods
        DESC: this method wraps all the steps used to
        store the stats .
        """
        try:
            pRate = self.qProfiler.get_profile_rate(self.previousProfile, self.playlist)
            sgmInfo = self.sStats.prepare_segment_stats(
                pRate, self.previousProfile, self.sgmUri, self.sgmDur)

            logicInfo = self.sStats.prepare_logic_stats(
                self.switchP, self.selectP, self.estRate)

            pBufInfo = self.sStats.prepare_pBuffer_stats(
                self.minPlayBuf, self.maxPlayBuf, self.curPlayBuf)

            self.statsDic = self.sStats.prepare_all_stats(
                self.statsDic, self.dwlStats, self.idleDwlTime, sgmInfo,
                logicInfo, pBufInfo, self.sgmIndex)

            self.sStats.write_stats(self.statsDic, self.statsFile)

        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None

    def est_bw_wrapper(self, iterTime, mName):
        """
        INPUT: iteration time (:int), estimation rate method name (:string)
        OUTPUT: 
        USES: the given as argument method
        DESC: This method is used for running the estimation rate methods.
        """
        try:
            avgRate = self.sgmRate
            varRate = 0
            prevRate = self.sgmRate
            while not self.KillThreads:
                if self.ssLogic == 2:
                    [self.estRate, avgRate, varRate] = getattr(Logic(), mName)(self.sgmRate, avgRate, varRate)
                if self.ssLogic == 3:
                    avgRate = sum(self.sgmRates) / len(self.sgmRates)
                    # avgRate = self.sgmRates/len(self.sgmRates)
                    self.estRate = getattr(Logic(), mName)(avgRate, prevRate)
                    prevRate = avgRate
                    # print self.estRate
                    # print avgRate

                # print self.estRate
                # print "CALCULATED %s" % time.ctime(time.time())
                time.sleep(iterTime)
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None

    def init_buf_state_operations(self):
        """
        INPUT: self
        OUTPUT:-
        USES: Segments, QualityProfiles and Statistics Classes methods
        DESC: This method performs the segment downloading,
        quality profile switching/selecting and saving stats
        operations for the initial buffering state.
        """
        try:

            if (self.variancy):
                ''' first Segment logic'''
                self.currentProfile = self.logics.first_time_l1(self.playlist)

            while self.sgmIndex <= self.initBufSize:
                '''Get selected segment's duration'''
                self.sgmDur = self.sSegment.get_segment_duration(
                    self.sgmIndex, self.currentProfile, self.playlist)

                '''Check if playout buffer size will exceed max threshold '''
                self.pBuffer.check_max_buffer_for_init_buf_state(
                    self.curPlayBuf['pBuf_3'], self.maxPlayBuf, self.sgmDur)

                '''Get selected segment's uri'''
                self.sgmUri = self.sSegment.get_segment_uri(
                    self.sgmIndex, self.currentProfile, self.playlist)
                '''Get selected segment'''
                self.dwlStats = self.sSegment.get_segment(self.sgmUri)
                '''Get selected segment's rate'''
                self.sgmRate = self.dwlStats.get('dwlRate')
                '''Get selected segment's download time'''
                self.dwlTime = self.dwlStats.get('dwlTime')
                ''' hold current profile's id'''
                self.previousProfile = self.currentProfile
                ''' append segment's  dwl rate'''
                self.sgmRates.append(self.sgmRate)

                ''' Calculate idle time between consecutive downloads '''
                if self.sgmIndex == 1:
                    self.idleDwlTime = 0
                    self.prevDwlTime = self.dwlStats.get('endDwlTime')
                else:
                    self.idleDwlTime = self.dwlStats.get('startDwlTime') - self.prevDwlTime
                    self.prevDwlTime = self.dwlStats.get('endDwlTime')

                ''' calculate the current playout buffer stats '''
                self.pBuffer.calc_curr_playout_buffer_for_init_buf_state(
                    self.curPlayBuf, self.sgmDur)

                if self.variancy:
                    '''simplest logic: selects as next profile
                    the one fitting best (<=) to the download rate of each segment'''
                    if self.ssLogic == 1:
                        self.estRate = self.sgmRate
                        [self.switchP, self.selectP, self.currentProfile] = self.logics.logic_1(self.sgmRate, self.currentProfile, self.playlist)

                    '''logic from the 'DISTRIBUTED & ADAPTIVE HTTP STREAMING' paper'''
                    if self.ssLogic == 2:
                        if self.firstTimeStartThread:
                            iterationTime = 2
                            self.estRate = self.sgmRate
                            self.avgRate = self.sgmRate
                            thread.start_new_thread(self.est_bw_wrapper, (iterationTime, 'est_bw_l2'))
                            self.firstTimeStartThread = False

                        [self.switchP, self.selectP, self.currentProfile] = self.logics.logic_1(self.estRate, self.currentProfile, self.playlist)

                    '''logic from the  'Improving Fairness in QoS and QoE domains for Adaptive Video Streaming' paper'''
                    if self.ssLogic == 3:
                        if self.firstTimeStartThread:
                            self.estRate = self.sgmRate
                            iterationTime = 2
                            thread.start_new_thread(self.est_bw_wrapper, (iterationTime, 'est_bw_l3'))
                            self.firstTimeStartThread = False

                        self.linkCapacity = self.sgmRate
                        [self.switchP, self.selectP, self.currentProfile] = self.logics.logic_2(self.linkCapacity, self.estRate, self.currentProfile, self.playlist)

                    ''' Adobe's OSMF logic as the  'QDASH: A QoE-aware DASH system' indicates' paper indicates'''
                    if self.ssLogic == 4:
                        self.estRate = self.sgmRate
                        [self.switchP, self.selectP, self.currentProfile] = self.logics.logic_3(self.dwlTime, self.sgmDur, self.currentProfile, self.playlist)

                    ''' logic from the 'QDASH: A QoE-aware DASH system' indicates' paper'''
                    if self.ssLogic == 5:
                        self.estRate = self.sgmRate
                        [self.switchP, self.selectP, self.currentProfile] = self.logics.logic_4(self.sgmRate, self.sgmDur, self.curPlayBuf['pBuf_3'], self.currentProfile, self.playlist)

                self.write_stats_wrapper()
                self.info()
                self.sgmIndex += 1
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None

    def playout_state_operations(self):
        """
        INPUT: self
        OUTPUT: returns flag indicating end of stream (:boolean)
        USES: Segments, QualityProfiles and Statistics Classes  methods
        DESC: This method performs the segment downloading,
        quality profile switching/selecting and saving stats
        operations for the Playout state.
        """
        try:

            while True:
                '''select and use buffering logic'''
                if self.bLogic == 1:
                    bFlag = True

                if bFlag:
                    '''Get selected segment's duration'''
                    self.sgmDur = self.sSegment.get_segment_duration(
                        self.sgmIndex, self.currentProfile, self.playlist)
                    '''Get selected segment's uri'''
                    self.sgmUri = self.sSegment.get_segment_uri(
                        self.sgmIndex, self.currentProfile, self.playlist)
                    '''Get selected segment'''
                    self.dwlStats = self.sSegment.get_segment(self.sgmUri)
                    '''Get selected segment's rate'''
                    self.sgmRate = self.dwlStats.get('dwlRate')
                    '''Get selected segment's download time'''
                    self.dwlTime = self.dwlStats.get('dwlTime')
                    ''' hold current profile's id'''
                    self.previousProfile = self.currentProfile
                    ''' append segment's  dwl rate'''
                    self.sgmRates.append(self.sgmRate)

                    ''' Calculate idle time between consecutive downloads '''
                    self.idleDwlTime = self.dwlStats.get('startDwlTime') - self.prevDwlTime
                    self.prevDwlTime = self.dwlStats.get('endDwlTime')

                    ''' calculate playout buffer stats '''
                    self.curPlayBuf = self.pBuffer.calc_curr_playout_buffer(
                        self.curPlayBuf, self.sgmDur, self.idleDwlTime, self.dwlStats['dwlTime'])

                    if self.variancy:
                        '''simplest logic: selects as next profile
                        the one fitting best (<=) to the download rate of each segment'''
                        if self.ssLogic == 1:
                            self.estRate = self.sgmRate
                            [self.switchP, self.selectP, self.currentProfile] = self.logics.logic_1(self.sgmRate, self.currentProfile, self.playlist)

                        '''logic from the 'DISTRIBUTED & ADAPTIVE HTTP STREAMING' paper'''
                        if self.ssLogic == 2:
                            [self.switchP, self.selectP, self.currentProfile] = self.logics.logic_1(self.estRate, self.currentProfile, self.playlist)

                        '''logic from the 'Improving Fairness in QoS and QoE domains for Adaptive Video Streaming' paper'''
                        if self.ssLogic == 3:
                            self.linkCapacity = self.sgmRate
                            [self.switchP, self.selectP, self.currentProfile] = self.logics.logic_2(self.linkCapacity, self.estRate, self.currentProfile, self.playlist)

                        ''' Adobe's OSMF logic as the  'QDASH: A QoE-aware DASH system' indicates' paper indicates'''
                        if self.ssLogic == 4:
                            self.estRate = self.sgmRate
                            [self.switchP, self.selectP, self.currentProfile] = self.logics.logic_3(self.dwlTime, self.sgmDur, self.currentProfile, self.playlist)

                        ''' logic from the 'QDASH: A QoE-aware DASH system' indicates' paper'''
                        if self.ssLogic == 5:
                            self.estRate = self.sgmRate
                            [self.switchP, self.selectP, self.currentProfile] = self.logics.logic_4(self.sgmRate, self.sgmDur, self.curPlayBuf['pBuf_3'], self.currentProfile, self.playlist)

                    self.write_stats_wrapper()
                    self.info()
                    self.sgmIndex += 1
                else:
                    pass

                '''break loop when no more segments exist)'''
                if self.sgmIndex > self.qProfiler.get_profile_segment_count(self.currentProfile, self.playlist):
                    break
        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None

    def manager(self, manifest, maxPlayBuf, minPlayBuf, initPBuf, ssLogic, bLogic):
        """
        INPUT: A manifest's uri (:string), max playout buffer (:int)
        min playout buffer (:int), initial buffering segments (:int),
        switch/select logic id (:int) and buffering logic (:int)
        OUTPUT: -
        USES:-
        DESC: This method manages all the DASH operations
        """
        try:

            ''' Manifest URI '''
            self.manifest = manifest
            ''' Initial Buffering Size (segments) '''
            self.initBufSize = initPBuf
            ''' Max playout Buffer Size '''
            self.maxPlayBuf = maxPlayBuf
            ''' Min playout Buffer Size '''
            self.minPlayBuf = minPlayBuf
            ''' Select/Switch Logic to use '''
            self.ssLogic = ssLogic
            ''' buffering Logic to use '''
            self.bLogic = bLogic

            ''' 
            current playout Buffer Size
            pBuf_1 is size at the start of segment download
            pBuf_2 is size at the end of the segment download
            not including the  downloaded segment
            pBuf_3 is size at the end of the segment download
            including the  downloaded segment
            '''
            self.curPlayBuf = OrderedDict(
                {"pBuf_1": 0.0, "pBuf_2": 0.0, "pBuf_3": 0.0})

            ''' statistics dictionary '''
            self.statsDic = OrderedDict()
            ''' Previous profile '''
            self.previousProfile = 1
            ''' Current profile '''
            self.currentProfile = 1
            ''' segment index '''
            self.sgmIndex = 1
            ''' Flag for ending Threads upon main program's end '''
            self.KillThreads = False
            ''' Flag for running Thread only once '''
            self.firstTimeStartThread = True
            ''' Bottleneck link's capacity '''
            self.linkCapacity = 10000000
            ''' hold sgm dwl rates'''
            self.sgmRates = []

            ''' JSON STATS FILE '''
            timestamp = round(time.time(), 6)
            self.statsFile = 'client_' + '%f' % timestamp + '_segmentStats.json'

            ''' Classes to use '''
            self.qProfiler = QualityProfiles()
            self.sSegment = Segments()
            self.sStats = Statistics()
            self.pBuffer = PlayoutBuffer()
            self.logics = Logic()

            ''' Check if max min values '''
            if self.maxPlayBuf < self.minPlayBuf:
                sys.exit("Min playout buffer value " +
                         "is larger than the max one. Please correct them!!!")

            ''' if min Playout buffer value is negative '''
            if self.minPlayBuf < 0:
                sys.exit("Min playout buffer value" +
                         "is negative. Please change it to a value >=0 !!!")

            ''' check if initial buffer segments is below 1 '''
            if self.initBufSize < 1:
                sys.exit("initial Buffer segments " +
                         "are below 1. Please change it to a value >=1 !!!")

            ''' MANIFEST PARSE OPERATIONS '''
            [self.playlist, self.variancy] = self.parse_operation()

            ''' CHECK IF WE HAVE A PLAYLIST '''
            if self.playlist is not None:
                print
                "We have a playlist. Let's play!!!"

                ''' Initialise the stats file '''
                with open(self.statsFile, 'w') as f:
                    json.dump({}, f)

                ''' START INITIAL BUFFERING '''
                self.init_buf_state_operations()

                ''' START PLAYOUT STATE '''
                self.playout_state_operations()

            else:
                sys.exit("Playlist of type None. Check parsing operation")

        except Exception as ex:
            print(sys._getframe().f_code.co_name + ": " + str(ex))
            return None
        finally:
            self.KillThreads = True

##testUri = "http://10.27.27.22/hls/bigBunny/segments/628kbps/628k.m3u8"
##testUri = "http://10.27.27.22/mpegDash/Sintel/segments/master.mpd"
##testUri = "http://10.27.27.22:8081/hls/bigBunny/segments/master-8081.m3u8"

##DashClientEmulator().manager(testUri, 30, 20, 3, 3, 1)
