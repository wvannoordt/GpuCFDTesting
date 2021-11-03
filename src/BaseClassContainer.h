#ifndef BASE_CLASS_CONTAINER_H
#define BASE_CLASS_CONTAINER_H
#include <vector>
#include <type_traits>
template <class BaseType> class BaseClassContainer
{
    public:
        
        BaseClassContainer(void){}
        
        template <typename T> void assert_all_copy_constructible(T t)
        {
            static_assert(std::is_trivially_copy_constructible<T>::value, "Argument type not trivially copy-constructible!");
        }
        template <typename T, typename... Ts> void assert_all_copy_constructible(T t, Ts... ts)
        {
            static_assert(std::is_trivially_copy_constructible<T>::value, "Argument type not trivially copy-constructible!");
            assert_all_copy_constructible(ts...);
        } 
        template <class NewObjType, typename... ts> NewObjType* Add(ts... params)
        {
            assert_all_copy_constructible(params...);
            static_assert(std::is_base_of<BaseType, NewObjType>::value, "Template type does not inherit from base type!");
            NewObjType* newObj = new NewObjType(params...);
            items.push_back(newObj);
            OnAdd(newObj);
            return newObj;
        }
        
        virtual void OnAdd(BaseType* newItem) {}
        
        ~BaseClassContainer(void)
        {
            for (auto& p:items)
            {
                delete p;
            }
        }
        
        template <class SearchType> std::vector<SearchType*> GetItemsByType(void)
        {
            std::vector<SearchType*> output;
            static_assert(std::is_base_of<BaseType, SearchType>::value, "Template type does not inherit from base type!");
            for (auto t:items)
            {
                SearchType* ptr = dynamic_cast<SearchType*>(t);
                if (ptr != NULL)
                {
                    output.push_back(ptr);
                }
            }
            return output;
        }

        BaseType* operator [] (int i) {return items[i];}
        
        typename std::vector<BaseType*>::iterator begin() noexcept
        {
            return items.begin();
        }

        typename std::vector<BaseType*>::const_iterator begin() const noexcept
        {
            return items.begin();
        }
        
        typename std::vector<BaseType*>::iterator end() noexcept
        {
            return items.end();
        }

        typename std::vector<BaseType*>::const_iterator end() const noexcept
        {
            return items.end();
        }
        
    protected:
        std::vector<BaseType*> items;
};

#endif