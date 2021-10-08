//
//  Service.swift
//  Weather App
//
//  Created by dato che on 2/5/21.
//

import Foundation

class Service {
    
    private let apiKey = "aaa274812b9713df7f2a9c4c39e635c8"
    private var components = URLComponents()
    
    init(for infoType: InfoType){
        components.scheme = "https"
        components.host = "api.openweathermap.org"
        components.path = "/data/2.5/" + infoType.rawValue
    }
    
    func loadWeather(for city: String, completion: @escaping (Result<ServiceResponse, Error>) -> ()){
        let parameters = [
            "q": city,
            "appid": apiKey.description
        ]
        loadWeather(for: parameters, completion: completion)
    }
    
    func loadWeather(for latitude: Double, longtitude: Double, completion: @escaping (Result<ServiceResponse, Error>) -> ()){
        let parameters = [
            "lat": latitude.description,
            "lon": longtitude.description,
            "appid": apiKey.description
        ]
        loadWeather(for: parameters, completion: completion)
    }
    
    private func loadWeather(for parameters: [String : String], completion: @escaping (Result<ServiceResponse, Error>) -> ()){
        components.queryItems = parameters.map{ key, value in
            return URLQueryItem(name: key, value: value)
        }
        
        if let url = components.url{
            let request = URLRequest(url: url)
            
            let task = URLSession.shared.dataTask(with: request, completionHandler: {data, response, error in
                if let error = error{
                    completion(.failure(error))
                    return
                }
                
                if let data = data{
                    let decoder = JSONDecoder()
                    do{
                        let result = try decoder.decode(ServiceResponse.self, from: data)
                        completion(.success(result))
                    } catch {
                        completion(.failure(ServiceError.noData))
                    }
                }
            })
            task.resume()
        }else{
            completion(.failure(ServiceError.invalidParameters))
        }
    
    }
    
    func loadForecast(for city: String, completion: @escaping (Result<FiveDayForecast, Error>) -> ()){
        let parameters = [
            "q": city,
            "appid": apiKey.description
        ]
        components.queryItems = parameters.map{ key, value in
            return URLQueryItem(name: key, value: value)
        }
        if let url = components.url{
            let request = URLRequest(url: url)
            
            let task = URLSession.shared.dataTask(with: request, completionHandler: {data, response, error in
                if let error = error{
                    completion(.failure(error))
                    return
                }
                
                if let data = data{
                    let decoder = JSONDecoder()
                    do{
                        let result = try decoder.decode(FiveDayForecast.self, from: data)
                        completion(.success(result))
                    } catch {
                        completion(.failure(ServiceError.noData))
                    }
                }
            })
            task.resume()
        }else{
            completion(.failure(ServiceError.invalidParameters))
        }
    }
}

enum ServiceError: Error{
    case noData
    case invalidParameters
}

enum InfoType: String{
    case today = "weather"
    case prolonged = "forecast"
}
